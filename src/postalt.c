/* biscuit postalt -- prefer primary chromosomes over unplaced/alt contigs.
 *
 * Reads an aligned SAM/BAM. For each primary record, if its own contig is an
 * unplaced/unlocalized/alt contig (--unplaced mode: names not matching the
 * primary-contig pattern) and its XA tag holds an equal/near-score hit on a
 * PRIMARY chromosome, the primary-chromosome hit is promoted to be the record
 * (the old placement is demoted into XA). MAPQ is then recomputed over distinct
 * PRIMARY loci only: unique-on-primary -> a confident value; >=2 primary loci
 * stay as they were (genuine ambiguity remains masked). Unplaced-only reads and
 * multi-primary reads are never boosted; MAPQ is never lowered.
 *
 * This is the unplaced-contig analogue of bwa-postalt, driven purely by XA; it
 * leaves `biscuit align` untouched. XA carries NM but no per-hit score, so the
 * near-score test uses NM (probes are fixed-length), and the restored MAPQ is a
 * fixed value (a clean unique hit is a true 60 anyway; downstream masking only
 * tests MAPQ thresholds).
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2026 Wanding.Zhou@pennmedicine.upenn.edu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include <unistd.h>
#include <getopt.h>
#include <regex.h>
#include "wzmisc.h"
#include "sam.h"
#include "hts.h"
#include "bamfilter.h"
#include "bisc_utils.h"

#define DEFAULT_PRIMARY_RE "^chr([0-9]{1,2}|X|Y|M|MT)$"

typedef struct {
  regex_t   pri_re;        // primary-contig pattern
  int8_t   *cls;           // per-tid: 1 primary, 0 secondary (lazy, sized n_targets)
  int       n_cls;
  int       nm_delta;      // a primary hit within bestNM+delta is "good"
  int       mapq_unique;   // MAPQ to set when unique-on-primary
  int       unplaced;      // mode flag (only mode implemented)
  // stats
  long n_reads, n_lift, n_restore, n_multi, n_pure_unplaced;
} postalt_conf_t;

/* one candidate alignment (record itself or an XA entry) */
typedef struct {
  int   tid;               // -1 if the contig is not in the header
  hts_pos_t pos;           // 0-based
  int   is_rev;
  char  cigar[64];
  int   nm;
  int   is_pri;            // on a primary chromosome
} hit_t;

static inline char dna_comp(char c) {
  switch (c) {
    case 'A': return 'T'; case 'T': return 'A';
    case 'C': return 'G'; case 'G': return 'C';
    case 'a': return 't'; case 't': return 'a';
    case 'c': return 'g'; case 'g': return 'c';
    default:  return 'N';
  }
}

/* classify a contig name as primary (1) or not (0) */
static int name_is_primary(postalt_conf_t *c, const char *name) {
  return regexec(&c->pri_re, name, 0, NULL, 0) == 0;
}

/* lazily build the per-tid primary/secondary class array from the header */
static void ensure_cls(postalt_conf_t *c, bam_hdr_t *h) {
  if (c->cls) return;
  c->n_cls = h->n_targets;
  c->cls = (int8_t*) malloc(c->n_cls * sizeof(int8_t));
  int i;
  for (i = 0; i < c->n_cls; ++i)
    c->cls[i] = name_is_primary(c, h->target_name[i]) ? 1 : 0;
}

/* parse a CIGAR string into a uint32 array; returns n_cigar (0 on '*'/error) */
static int parse_cigar(const char *s, uint32_t *cig, int max_cig) {
  int n = 0; long len = 0;
  if (!s || s[0] == '*') return 0;
  for (; *s; ++s) {
    if (*s >= '0' && *s <= '9') { len = len*10 + (*s - '0'); }
    else {
      int op = bam_cigar_table[(int)*s];
      if (op < 0 || n >= max_cig) return 0;
      cig[n++] = bam_cigar_gen((uint32_t)len, op);
      len = 0;
    }
  }
  return n;
}

/* format a bam1_t's own CIGAR into a string buffer */
static void record_cigar_str(bam1_t *b, char *out, int cap) {
  uint32_t *c = bam_get_cigar(b); int i, k = 0;
  if (b->core.n_cigar == 0) { out[0]='*'; out[1]=0; return; }
  for (i = 0; i < (int)b->core.n_cigar && k < cap-12; ++i)
    k += snprintf(out+k, cap-k, "%u%c", bam_cigar_oplen(c[i]),
                  "MIDNSHP=X"[bam_cigar_op(c[i])]);
}

/* rebuild a record at a new placement (handles strand flip / cigar change) */
static void rebuild_at(bam1_t *b, bam_hdr_t *h, int tid, hts_pos_t pos,
                       int is_rev, uint32_t *cig, int ncig, uint8_t mapq) {
  int lq = b->core.l_qseq;
  int cur_rev = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
  /* capture aux blob, seq (ascii, reoriented), qual (reoriented) */
  uint8_t *aux = bam_get_aux(b);
  int l_aux = b->l_data - (int)(aux - b->data);
  uint8_t *auxcpy = (uint8_t*) malloc(l_aux > 0 ? l_aux : 1);
  memcpy(auxcpy, aux, l_aux);

  char *seq  = (char*) malloc(lq + 1);
  char *qual = (char*) malloc(lq > 0 ? lq : 1);
  uint8_t *s = bam_get_seq(b), *q = bam_get_qual(b);
  int no_qual = (lq > 0 && q[0] == 0xff);
  int i;
  if (is_rev == cur_rev) {
    for (i = 0; i < lq; ++i) { seq[i] = seq_nt16_str[bam_seqi(s,i)]; qual[i] = q[i]; }
  } else {
    for (i = 0; i < lq; ++i) {
      seq[lq-1-i]  = dna_comp(seq_nt16_str[bam_seqi(s,i)]);
      qual[lq-1-i] = q[i];
    }
  }
  seq[lq] = 0;

  uint16_t flag = b->core.flag;
  if (is_rev) flag |= BAM_FREVERSE; else flag &= ~BAM_FREVERSE;

  char *qname = strdup(bam_get_qname(b));
  bam_set1(b, strlen(qname), qname, flag, tid, pos, mapq, ncig, cig,
           b->core.mtid, b->core.mpos, b->core.isize,
           lq, seq, no_qual ? NULL : qual, l_aux);
  memcpy(bam_get_aux(b), auxcpy, l_aux);   // reserved region filled with old aux

  free(auxcpy); free(seq); free(qual); free(qname);
}

/* build an XA string from all hits except index `skip` (the promoted one).
   hits[0] is the record's old placement, so it is naturally demoted here.
   NM-only, biscuit format. */
static void build_xa(kstring_t *xa, hit_t *hits, int nhit, int skip,
                     bam_hdr_t *h) {
  xa->l = 0;
  int i;
  for (i = 0; i < nhit; ++i) {
    if (i == skip) continue;
    if (xa->l) kputc(';', xa);
    ksprintf(xa, "%s,%c%ld,%s,%d", h->target_name[hits[i].tid],
             "+-"[hits[i].is_rev], (long)(hits[i].pos + 1), hits[i].cigar, hits[i].nm);
  }
}

static int postalt_func(bam1_t *b, samFile *out, bam_hdr_t *h, void *data) {
  postalt_conf_t *c = (postalt_conf_t*) data;
  ensure_cls(c, h);

  if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
    if (out && sam_write1(out, h, b) < 0) wzfatal("Cannot write bam.\n");
    return 0;
  }
  c->n_reads++;

  /* gather hits: the record itself + XA entries */
  hit_t hits[256]; int nhit = 0;
  uint8_t *nm_tag = bam_aux_get(b, "NM");
  hits[0].tid = b->core.tid;
  hits[0].pos = b->core.pos;
  hits[0].is_rev = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
  record_cigar_str(b, hits[0].cigar, sizeof(hits[0].cigar));
  hits[0].nm = nm_tag ? bam_aux2i(nm_tag) : 0;
  hits[0].is_pri = (b->core.tid < c->n_cls) ? c->cls[b->core.tid] : 0;
  nhit = 1;

  uint8_t *xa_tag = bam_aux_get(b, "XA");
  char *xa_str = xa_tag ? bam_aux2Z(xa_tag) : NULL;
  char *xa_dup = xa_str ? strdup(xa_str) : NULL;
  if (xa_dup) {
    char *save1, *e;
    for (e = strtok_r(xa_dup, ";", &save1); e && nhit < 256;
         e = strtok_r(NULL, ";", &save1)) {
      /* entry: name,±pos,CIGAR,NM */
      char *save2, *tok;
      char *name = strtok_r(e, ",", &save2);
      char *spos = strtok_r(NULL, ",", &save2);
      char *cig  = strtok_r(NULL, ",", &save2);
      char *nm   = strtok_r(NULL, ",", &save2);
      if (!name || !spos || !cig || !nm) continue;
      int tid = bam_name2id(h, name);
      if (tid < 0) continue;
      hit_t *ht = &hits[nhit];
      ht->tid = tid;
      ht->is_rev = (spos[0] == '-');
      ht->pos = (hts_pos_t) strtol(spos + 1, NULL, 10) - 1; // 1-based -> 0-based
      snprintf(ht->cigar, sizeof(ht->cigar), "%s", cig);
      ht->nm = atoi(nm);
      ht->is_pri = (tid < c->n_cls) ? c->cls[tid] : 0;
      (void)tok;
      nhit++;
    }
  }

  /* best NM and "good" primary hits (within bestNM + delta) */
  int bestNM = hits[0].nm, i;
  for (i = 1; i < nhit; ++i) if (hits[i].nm < bestNM) bestNM = hits[i].nm;

  /* best near-equal primary hit -> the promotion target (uses nm_delta) */
  int best_pri = -1;
  for (i = 0; i < nhit; ++i) {
    if (!hits[i].is_pri || hits[i].nm > bestNM + c->nm_delta) continue;
    if (best_pri < 0 || hits[i].nm < hits[best_pri].nm) best_pri = i;
  }

  /* mask decision: count DISTINCT PRIMARY loci at ANY score (biscuit already
     filtered XA to near-best via its drop ratio, so every listed primary is a
     real competitor), plus the number of unplaced hits. We only restore MAPQ
     with positive evidence: exactly one primary locus AND >=1 unplaced hit that
     explains biscuit's low MAPQ. An empty XA (n_unplaced==0) means the low MAPQ
     came from something biscuit did not list -- leave it masked. */
  int n_pri_loci = 0, n_unplaced = 0;
  int loc_tid[256]; long loc_pos[256];
  int win = b->core.l_qseq > 0 ? b->core.l_qseq : 50;
  for (i = 0; i < nhit; ++i) {
    if (!hits[i].is_pri) { n_unplaced++; continue; }
    int j, dup = 0;
    for (j = 0; j < n_pri_loci; ++j)
      if (loc_tid[j] == hits[i].tid && labs((long)hits[i].pos - loc_pos[j]) <= win) { dup = 1; break; }
    if (!dup && n_pri_loci < 256) { loc_tid[n_pri_loci] = hits[i].tid; loc_pos[n_pri_loci] = hits[i].pos; n_pri_loci++; }
  }

  int rec_pri = hits[0].is_pri;
  int did_lift = 0;

  /* (1) promote: record on a secondary contig, a good primary exists */
  if (!rec_pri && best_pri > 0) {
    hit_t *p = &hits[best_pri];
    kstring_t newxa = {0,0,0};
    build_xa(&newxa, hits, nhit, best_pri, h);   // includes hits[0] (old placement)

    uint32_t cig[64];
    int ncig = parse_cigar(p->cigar, cig, 64);
    uint8_t new_mapq = b->core.qual;   // mapq decided below; keep for now
    if (ncig == 0) {
      /* unparseable cigar -> leave record untouched, just write */
      free(newxa.s);
    } else if (p->is_rev == hits[0].is_rev &&
               ncig == (int)b->core.n_cigar) {
      /* fast in-place: same strand, same cigar length */
      b->core.tid = p->tid;
      b->core.pos = p->pos;
      uint32_t *rc = bam_get_cigar(b);
      int k; for (k = 0; k < ncig; ++k) rc[k] = cig[k];
      bam_aux_update_int(b, "NM", p->nm);
      if (newxa.l) bam_aux_update_str(b, "XA", newxa.l + 1, newxa.s);
      else { uint8_t *xa = bam_aux_get(b, "XA"); if (xa) bam_aux_del(b, xa); }
      did_lift = 1;
      free(newxa.s);
    } else {
      /* full rebuild (strand flip or cigar change) */
      rebuild_at(b, h, p->tid, p->pos, p->is_rev, cig, ncig, new_mapq);
      bam_aux_update_int(b, "NM", p->nm);
      if (newxa.l) bam_aux_update_str(b, "XA", newxa.l + 1, newxa.s);
      else { uint8_t *xa = bam_aux_get(b, "XA"); if (xa) bam_aux_del(b, xa); }
      did_lift = 1;
      free(newxa.s);
    }
    if (did_lift) c->n_lift++;
  }

  /* (2) restore MAPQ only with positive evidence the ambiguity is unplaced:
     the (possibly lifted) record sits on the sole primary locus AND >=1 unplaced
     hit explains biscuit's low MAPQ. Never restore on an empty XA or when any
     other primary locus competes. */
  int final_on_primary = rec_pri || did_lift;
  if (final_on_primary && n_pri_loci == 1 && n_unplaced >= 1) {
    if (b->core.qual < c->mapq_unique) { b->core.qual = c->mapq_unique; c->n_restore++; }
  } else if (n_pri_loci >= 2) {
    c->n_multi++;                       // genuine multi-primary: leave masked
  } else {
    c->n_pure_unplaced++;               // no unplaced evidence / unplaced-only: leave
  }

  free(xa_dup);
  if (out && sam_write1(out, h, b) < 0) wzfatal("Cannot write bam.\n");
  return 0;
}

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: biscuit postalt [options] <in.bam> [out.bam]\n");
  fprintf(stderr, "Prefer primary chromosomes over unplaced/alt contigs using XA.\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    --unplaced       enable unplaced/alt-contig demotion (the only mode)\n");
  fprintf(stderr, "    -p STR           primary-contig regex [%s]\n", DEFAULT_PRIMARY_RE);
  fprintf(stderr, "    -d INT           NM near-score delta (a primary hit within\n");
  fprintf(stderr, "                     bestNM+delta is a promotion candidate) [1]\n");
  fprintf(stderr, "    -q INT           MAPQ to set when unique on the primary assembly [60]\n");
  fprintf(stderr, "    -g STR           process only this region (needs an index)\n");
  fprintf(stderr, "    -h               this help\n\n");
  fprintf(stderr, "in.bam/out.bam may be '-' for stdin/stdout; input may be SAM or BAM.\n\n");
  return 1;
}

int main_postalt(int argc, char *argv[]) {
  postalt_conf_t c = {0};
  c.nm_delta = 1; c.mapq_unique = 60; c.unplaced = 0;
  char *pri_pat = DEFAULT_PRIMARY_RE;
  char *reg = 0;

  kstring_t cl = generate_command_line_string(argc, argv);

  static struct option lopt[] = {
    {"unplaced", no_argument, 0, 1000},
    {0,0,0,0}
  };
  int c0;
  while ((c0 = getopt_long(argc, argv, ":p:d:q:g:h", lopt, NULL)) >= 0) {
    switch (c0) {
      case 1000: c.unplaced = 1; break;
      case 'p': pri_pat = optarg; break;
      case 'd': c.nm_delta = atoi(optarg); break;
      case 'q': c.mapq_unique = atoi(optarg); break;
      case 'g': reg = optarg; break;
      case 'h': return usage();
      case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
      case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
      default: return usage();
    }
  }
  if (!c.unplaced) { usage(); wzfatal("This build only supports --unplaced mode.\n"); }

  char *infn  = optind < argc ? argv[optind++] : NULL;
  char *outfn = optind < argc ? argv[optind++] : "-";
  if (!infn) { usage(); wzfatal("Please provide an input SAM/BAM (or '-').\n"); }

  if (regcomp(&c.pri_re, pri_pat, REG_EXTENDED | REG_NOSUB) != 0)
    wzfatal("Bad primary-contig regex: %s\n", pri_pat);

  bam_filter(infn, outfn, reg, &c, cl.s, postalt_func);

  fprintf(stderr, "[postalt] reads=%ld lifted=%ld mapq-restored=%ld "
          "multi-primary(masked)=%ld unplaced-only=%ld\n",
          c.n_reads, c.n_lift, c.n_restore, c.n_multi, c.n_pure_unplaced);

  regfree(&c.pri_re);
  free(c.cls);
  free(cl.s);
  return 0;
}
