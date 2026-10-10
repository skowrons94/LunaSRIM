"""Random number generator of TRIM (SRIM-2013).

TRIM uses L'Ecuyer's combined multiplicative congruential generator
(the same as Numerical Recipes' ``ran2`` but *without* the Bays-Durham
shuffle).  The integer division is performed with Visual Basic's ``CLng``
which rounds to the nearest (even) integer instead of truncating.  The two
seeds are taken from the "Random Number Seed" line of TRIM.IN (a zero seed is
replaced by 16381 and 12345 respectively).
"""

M1 = 2147483563
M2 = 2147483399


def _clng(x):
    # Visual Basic CLng: round half to even
    return int(round(x))


class TrimRandom(object):
    __slots__ = ('s1', 's2', 'count')

    def __init__(self, seed1=0, seed2=0):
        self.seed(seed1, seed2)

    def seed(self, seed1=0, seed2=0):
        self.s1 = int(seed1) if seed1 else 16381
        self.s2 = int(seed2) if seed2 else 12345
        self.count = 0

    def random(self):
        s1 = self.s1
        k = _clng(s1 / 53668.0)
        s1 = 40014 * (s1 - k * 53668) - k * 12211
        if s1 < 0:
            s1 += 2147483563
        s2 = self.s2
        k = _clng(s2 / 52774.0)
        s2 = 40692 * (s2 - k * 52774) - k * 3791
        if s2 < 0:
            s2 += 2147483399
        self.s1 = s1
        self.s2 = s2
        z = s1 - s2
        if z < 1:
            z += 2147483562
        self.count += 1
        return z / 2147483563.0

    __call__ = random

    def getstate(self):
        return (self.s1, self.s2, self.count)

    def setstate(self, st):
        self.s1, self.s2, self.count = st
