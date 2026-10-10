"""SRIM's compound dictionary (Data/Compound.dat): densities and
stoichiometries of ~350 common targets (gases, plastics, biological
materials, implantation compounds...)."""
import csv
import re

from . import data_files


class Compound(object):
    def __init__(self, category, title, name, density, elements, gas, mass_percent):
        self.category = category
        self.title = title            # description shown in SRIM's menu
        self.name = name
        self.density = density        # g/cm3
        self.elements = elements      # [(Z, atom fraction)], fractions normalised to 1
        self.gas = gas
        self.mass_percent = mass_percent

    def as_dict(self):
        return {'category': self.category, 'title': self.title, 'name': self.name,
                'density': self.density, 'gas': self.gas,
                'elements': [{'z': z, 'fraction': f} for z, f in self.elements]}


def _category(line):
    if '\xdf' not in line:            # section titles are framed with this character
        return None
    text = re.sub(r'[^A-Za-z0-9/ ,()&.-]', ' ', line[1:]).strip()
    text = re.sub(r'\s+', ' ', text)
    if not text or any(w in text for w in ('Density', 'Name', 'Stoich', 'g/cm3')):
        return None
    text = text.split(' from ')[0].split('(9')[0].strip()
    return text.title().replace('(Misc.)', '(misc.)')


_CACHE = {}


def load(directory=data_files.DATA_DIR):
    if directory in _CACHE:
        return _CACHE[directory]
    path = data_files._find('Compound.dat', directory)
    with open(path, 'rb') as fh:
        lines = fh.read().decode('latin-1').replace('\r\n', '\n').split('\n')
    atoms = data_files.load_atomdata(directory)
    out = []
    cat = 'Compounds'
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith('!'):
            c = _category(l)
            if c:
                cat = c
        elif l.startswith('*') and i + 1 < len(lines):
            title = re.split(r'\s{2,}', l[1:].strip())[0]
            gas = '(gas)' in l.lower()
            try:
                row = next(csv.reader([lines[i + 1]], skipinitialspace=True))
                name = row[0].strip()
                density = float(row[1])
                n = int(float(row[2]))
                vals = [float(x) for x in row[3:3 + 2 * n] if x.strip()]
                if len(vals) < 2 * n:
                    raise ValueError(name)
            except (ValueError, IndexError, StopIteration):
                i += 1
                continue
            mass_pc = name.startswith('%')
            els = []
            for k in range(n):
                z, amount = int(vals[2 * k]), vals[2 * k + 1]
                if mass_pc:
                    amount = amount / atoms[z].natural_weight
                els.append((z, amount))
            tot = sum(a for z, a in els) or 1.0
            els = [(z, a / tot) for z, a in els if a > 0]
            gas = gas or density < 0.01
            out.append(Compound(cat, title, name.lstrip('%'), density, els, gas, mass_pc))
            i += 1
        i += 1
    _CACHE[directory] = out
    return out
