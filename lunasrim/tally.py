"""Accumulation of TRIM's distributions (ionization, phonons, vacancies, ...).

The arrays use TRIM's own layout so that they can be compared element by
element with the program: 2D arrays are stored flat with the first (VB)
index running fastest, and all indices are TRIM's 1-based ones (atom 1 is
the ion, atoms 2..n+1 are the target elements, layers are 1..nlay).
"""


class VBArray(object):
    """flat VB array with given (lbound, count) per dimension"""
    __slots__ = ('name', 'dims', 'data', 'log')

    def __init__(self, name, dims, log=None):
        self.name = name
        self.dims = dims
        n = 1
        for lb, cnt in dims:
            n *= cnt
        self.data = [0] * n          # VB Empty -> 0
        self.log = log

    def index(self, i, j=None):
        (lb1, n1) = self.dims[0]
        if j is None:
            return i - lb1
        lb2 = self.dims[1][0]
        return (i - lb1) + n1 * (j - lb2)

    def get(self, i, j=None):
        return self.data[self.index(i, j)]

    def set(self, v, i, j=None):
        k = self.index(i, j)
        self.data[k] = v
        if self.log is not None:
            self.log.write('%s[%d] = %r\n' % (self.name, k, v))

    def add(self, v, i, j=None):
        """element = element + v"""
        k = self.index(i, j)
        r = self.data[k] + v
        self.data[k] = r
        if self.log is not None:
            self.log.write('%s[%d] = %r\n' % (self.name, k, r))


def vb_int(x):
    """VB Int(): floor"""
    import math
    return float(math.floor(x))


class TrimTally(object):
    """The distribution arrays of TRIM (names are TRIM's global addresses)."""

    def __init__(self, nel, nlay, log=None):
        self.log = log
        A = lambda name, *dims: VBArray(name, dims, log)
        n1 = nel + 1
        self.d95c = A('D4fd95c', (1, n1), (1, 100))     # vacancies/energy by atom & depth
        self.d964 = A('D4fd964', (1, 2), (1, 100))      # ionization: 1 ion, 2 recoils
        self.d968 = A('D4fd968', (1, n1))
        self.d96c = A('D4fd96c', (1, 100))
        self.d970 = A('D4fd970', (1, 2), (1, 100))      # phonons: 1 ion, 2 recoils
        self.d974 = A('D4fd974', (1, 100), (1, 100))
        self.d978 = A('D4fd978', (1, 100), (1, 100))
        self.d97c = A('D4fd97c', (1, 100), (1, 100))    # ionization depth x lateral
        self.d980 = A('D4fd980', (1, 100), (1, 100))    # phonons depth x lateral
        self.d984 = A('D4fd984', (1, 100), (1, 100))
        self.d988 = A('D4fd988', (1, 100), (1, 100))    # vacancies depth x lateral
        self.d98c = A('D4fd98c', (1, n1), (1, 100))
        self.d994 = A('D4fd994', (1, n1))
        self.d9a0 = A('D4fd9a0', (1, n1), (1, 100))     # vacancies by atom & depth
        self.d9ac = A('D4fd9ac', (1, 2), (1, 100))
        self.d9b0 = A('D4fd9b0', (1, 2), (1, 100))
        self.d9b8 = A('D4fd9b8', (1, n1), (1, 200))
        # scalar accumulators
        self.g244 = 0.0      # vacancies
        self.g61c = 0.0      # displacements
        self.g230 = 0.0      # sputtered atoms
        self.g23c = 0.0
        self.g1c4 = 0.0
        self.g1d8 = 0.0      # replacement collisions
        self.g280 = self.g288 = self.g290 = self.g298 = 0.0   # range moments
        self.g2c0 = self.g2c8 = self.g2d0 = self.g2d8 = 0.0   # lateral moments
        self.n_back = 0
        self.n_trans = 0
        # running statistics (0x45a200)
        self.g_a38 = self.g_a48 = self.g_a10 = self.g_a20 = self.g_9f8 = self.g_a08 = 0
