"""Readers for the SRIM-2013 data files (SCOEF03.DAT, SNUC03.dat).

The column layout of SCOEF03.DAT was recovered from the SRIM SR-module
loader routine.  Every row describes one element (Z = 1..92) and contains two
quoted strings followed by 93 numbers:

====== ========================================================================
cols   meaning
====== ========================================================================
1      Z
2      mass of most abundant isotope (integer, amu)
3      mass of most abundant isotope (amu)
4      natural (isotope averaged) atomic weight (amu)
5      density of the solid (g/cm3)
6      atomic density of the solid (atoms/cm3)
7      Fermi velocity of the solid (v0 units)
8      heat of sublimation (eV)
9      density of the gas (g/cm3)
10     atomic density of the gas (atoms/cm3)
11-20  proton stopping fit coefficients            (``P[1..10]``)
21-30  helium stopping fit coefficients            (``HE[1..10]``)
31-40  lithium stopping fit coefficients           (``LI[1..10]``)
41-53  heavy ion effective charge polynomial  (ion) (``ION_POLY[1..13]``)
54-66  heavy ion scaling polynomial        (target) (``TGT_POLY[1..13]``)
67     mean ionisation potential I (eV), solid
68-69  mean ionisation potential I (eV), gas / alternative
70-80  shell correction polynomial                  (``SHELL_C[1..11]``)
81-88  shell correction polynomial                  (``SHELL_D[1..8]``)
89-93  high energy heavy-ion effective charge parameters (``HI[1..5]``)
====== ========================================================================
"""
import os
import shlex

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def _find(name, directory):
    p = os.path.join(directory, name)
    if os.path.exists(p):
        return p
    for f in os.listdir(directory):
        if f.lower() == name.lower():
            return os.path.join(directory, f)
    raise FileNotFoundError(name)


class ElementData(object):
    """Per-element data from SCOEF03.DAT."""

    def __init__(self, tokens):
        self.symbol = tokens[0].strip()
        self.name = tokens[1].strip()
        v = [float(x) for x in tokens[2:]]
        col = [None] + v                              # 1-based columns
        self.z = int(col[1])
        self.mai_mass_int = col[2]
        self.mai_mass = col[3]
        self.natural_mass = col[4]
        self.density = col[5]
        self.atomic_density = col[6]
        self.fermi_velocity = col[7]
        self.heat_sublimation = col[8]
        self.gas_density = col[9]
        self.gas_atomic_density = col[10]
        self.p_coef = col[11:21]
        self.he_coef = col[21:31]
        self.li_coef = col[31:41]
        self.ion_poly = col[41:54]
        self.tgt_poly = col[54:67]
        self.i_solid = col[67]
        self.i_gas = col[68:70]
        self.shell_c = col[70:81]
        self.shell_d = col[81:89]
        self.hi = col[89:94]


_SCOEF = {}
_SNUC = {}


def load_scoef(directory=DATA_DIR):
    if directory in _SCOEF:
        return _SCOEF[directory]
    table = {}
    with open(_find('SCOEF03.DAT', directory), encoding='latin-1') as fh:
        for line in fh:
            t = shlex.split(line)
            if len(t) < 95:
                continue
            e = ElementData(t)
            table[e.z] = e
    _SCOEF[directory] = table
    return table


def load_snuc(directory=DATA_DIR):
    """Nuclear stopping fit coefficients for every (Z1, Z2) pair."""
    if directory in _SNUC:
        return _SNUC[directory]
    table = {}
    with open(_find('SNUC03.dat', directory), encoding='latin-1') as fh:
        for line in fh:
            t = line.split()
            if len(t) < 6:
                continue
            z1, z2 = int(t[0]), int(t[1])
            c = [float(x) for x in t[2:6]]
            table[(z1, z2)] = c
            table[(z2, z1)] = c
    _SNUC[directory] = table
    return table


def element(z, directory=DATA_DIR):
    return load_scoef(directory)[int(z)]


def symbol_to_z(symbol, directory=DATA_DIR):
    s = symbol.strip().lower()
    for z, e in load_scoef(directory).items():
        if e.symbol.lower() == s:
            return z
    raise KeyError(symbol)
