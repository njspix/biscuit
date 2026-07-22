import logging
import sys

logger = logging.getLogger(__name__)

class ReadPair:
    """Store data from paired-end reads to check for correct sequence alignment.

    Args:
        chrm (str): known chromosome for alignment (from read name)
        start (int): known start position for alignment (from read name)
        end (int): known end position for alignment (from read name)
    """
    def __init__(self, chrm, start, end):
        self.chrm = chrm
        self.start = start
        self.end = end
        self.r1 = None
        self.r2 = None

    def __str__(self):
        return f'{self.chrm}:{self.start}-{self.end} -- R1: {self.r1} || R2: {self.r2}'

    def validate(self, chrm, start, end):
        """Validate known alignment matches initialized values. Must be run prior to adding new
        read to initialized ReadPair.

        Args:
            chrm (str): known chromosome for alignment (from read name)
            start (int): known start position for alignment (from read name)
            end (int): known end position for alignment (from read name)

        Returns:
            bool
        """
        return self.chrm == chrm and self.start == start and self.end == end

    def _create_dict(self, chrm, qloc, tlen):
        """Helper function for storing information from aligned read.

        Args:
            chrm (str): aligned chromosome
            qloc (str): aligned location
            tlen (str): aligned template length

        Returns:
            dict
        """
        return {
            'chrm': chrm,
            'qloc': int(qloc),
            'tlen': int(tlen),
        }

    def add_r1(self, chrm, qloc, tlen):
        """Add aligned read 1 (based on read name) to initialized ReadPair. Modifies existing
        ReadPair.

        Args:
            chrm (str): aligned chromosome
            qloc (str): aligned location (converted to int)
            tlen (str): aligned template length (converted to int)

        Returns:
            None
        """
        self.r1 = self._create_dict(chrm, qloc, tlen)

    def add_r2(self, chrm, qloc, tlen):
        """Add aligned read 2 (based on read name) to initialized ReadPair. Modifies existing
        ReadPair.

        Args:
            chrm (str): aligned chromosome
            qloc (str): aligned location (converted to int)
            tlen (str): aligned template length (converted to int)

        Returns:
            None
        """
        self.r2 = self._create_dict(chrm, qloc, tlen)

    def get_actual_region(self):
        """Get known region read pair came from. Region is defined as the start of the leftmost
        read to the end of the rightmost read.

        Args:
            None

        Returns:
            tuple
        """
        return (self.start, self.end)

    def get_mapped_region(self):
        """Get mapped region from added reads. Region is defined as the start of the leftmost read
        to the end of the rightmost read. Uses the template length (TLEN) value to determine region
        width.

        Args:
            None

        Returns:
            tuple
        """
        if self.r1 is None or self.r2 is None:
            logger.error('Read 1 or Read 2 is missing from ReadPair')
            sys.exit(1)

        # -1 is to account for 0-based vs 1-based regions
        if self.r1['tlen'] >= 0 and self.r2['tlen'] >= 0:
            # On the off chance that inputs are from BISCUIT v1.9.0 or earlier
            # Leftmost read needs to be figured out via comparison, not sign check of TLEN
            left = self.r1['qloc'] if self.r1['qloc'] < self.r2['qloc'] else self.r2['qloc']
            return (left, left + self.r1['tlen'] - 1)
        elif self.r1['tlen'] >= 0 and self.r2['tlen'] < 0:
            return (self.r1['qloc'], self.r1['qloc'] + self.r1['tlen'] - 1)
        else:
            # Technically, this encompasses two negative TLEN values, but I don't think this
            # can occur, so leave as is until it becomes a problem
            return (self.r2['qloc'], self.r2['qloc'] + self.r2['tlen'] - 1)

def name_to_pieces(name):
    """Turn Sherman read name into its constituent parts: read number, true chromosome, true start,
    true end, and read in pair.

    Args:
        name (str): Sherman read name

    Returns:
        tuple
    """
    pieces = name.split('_')

    # Initial parsing
    n = pieces[0] # read number
    pos = pieces[1] # chr:start-end
    read = pieces[2] # read in pair

    # True location parsing
    x_pos = pos.split(':')

    chr = x_pos[0]
    beg = int(x_pos[1].split('-')[0])
    end = int(x_pos[1].split('-')[1])

    return (n, read, chr, beg, end)

def read_file(fname):
    """Read SAM file and find number of misalignments.

    Args:
        fname (str): SAM file name

    Returns:
        tuple
    """
    data = {}
    with open(fname, 'r') as fh:
        for line in fh.readlines():
            # Skip header lines
            if line.startswith('@'):
                continue

            # Parse read entry
            l = line.strip().split('\t')
            n, read, c_actual, s_actual, e_actual = name_to_pieces(l[0])
            chrm = l[2]
            qloc = l[3]
            tlen = l[8]

            # If read in data already, validate names match; otherwise, create new entry
            if n in data:
                if not data[n].validate(c_actual, s_actual, e_actual):
                    logger.error(f'True locations do not match for read {n}')
                    sys.exit(1)
            else:
                data[n] = ReadPair(c_actual, s_actual, e_actual)

            # Add data
            if read == 'R1':
                data[n].add_r1(chrm, qloc, tlen)
            else:
                data[n].add_r2(chrm, qloc, tlen)

    # Compare true vs mapped locations
    n_diff = 0
    for n, pair in data.items():
        r_actual = pair.get_actual_region()
        r_mapped = pair.get_mapped_region()

        if not r_actual[0] == r_mapped[0] and not r_actual[1] == r_mapped[1]:
            logger.debug(f'{n}\t{pair}')
            n_diff += 1

    return n_diff, len(data)

def main(f1, f2):
    """Check alignments of two inputs files.

    Args:
        old (str): Path of old file
        new (str): Path of new file

    Returns:
        None
    """
    wrong_1, n_1 = read_file(f1)
    wrong_2, n_2 = read_file(f2)

    logger.info(f'OLD --- Incorrect mapping rate: {wrong_1/n_1:.3f} ({wrong_1} / {n_1}) --- File: {f1}')
    logger.info(f'NEW --- Incorrect mapping rate: {wrong_2/n_2:.3f} ({wrong_2} / {n_2}) --- File: {f2}')

    return None

if __name__ == '__main__':
    main()
