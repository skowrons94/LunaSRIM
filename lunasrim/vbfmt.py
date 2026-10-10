"""Visual Basic compatible number formatting (Str$, CStr, Format$).

TRIM writes its output files with Visual Basic 5 run-time functions; the
helpers below reproduce the conversions used there:

* ``vb_str``     ``Str$(x)``: leading blank for non-negative numbers, no
  leading zero for |x| < 1, at most 15 significant digits.
* ``vb_cstr``    ``CStr(x)`` (same without the sign blank).
* ``vb_format``  ``Format$(x, fmt)`` for numeric formats made of ``0``,
  ``#``, ``.``, ``E+``/``E-`` and literal characters, with optional
  positive;negative sections.
* ``vb_print_num`` how ``Print #`` writes a number (Str$ plus a blank).
"""
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN


def _digits15(x):
    """decimal representation of x with 15 significant digits"""
    return Decimal('%.14e' % x)


def _trim(s):
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s


def vb_cstr(x, single=False):
    """CStr() of a Double/Single/Integer (single: x is a VB Single, 7 digits)"""
    if isinstance(x, bool):
        return 'True' if x else 'False'
    if isinstance(x, int):
        return str(x)
    if x == 0:
        return '0'
    d = Decimal('%.6e' % x) if single else _digits15(x)
    sign = '-' if d < 0 else ''
    d = abs(d)
    exp = d.adjusted()
    if exp >= (7 if single else 15) or exp < -4:
        m = d.scaleb(-exp)
        ms = _trim(format(m, 'f'))
        return '%s%sE%s%02d' % (sign, ms, '+' if exp >= 0 else '-', abs(exp))
    s = _trim(format(d, 'f'))
    if s.startswith('0.'):
        s = s[1:]
    return sign + s


def vb_str(x, single=False):
    """Str$(): like CStr with a leading blank for positive numbers"""
    s = vb_cstr(x, single)
    if not s.startswith('-'):
        s = ' ' + s
    return s


def vb_print_num(x):
    """item text of a number written with Print #"""
    return vb_str(x) + ' '


def _round(d, ndec):
    q = Decimal(1).scaleb(-ndec)
    return d.quantize(q, rounding=ROUND_HALF_UP)


def _split_sections(fmt):
    secs = []
    cur = ''
    inq = False
    for ch in fmt:
        if ch == '"':
            inq = not inq
        if ch == ';' and not inq:
            secs.append(cur)
            cur = ''
        else:
            cur += ch
    secs.append(cur)
    return secs


def _parse(sec):
    """split a format section into prefix literal, numeric spec, suffix"""
    i = 0
    n = len(sec)
    out = []
    # tokens: ('lit', text) or ('num', spec)
    num = None
    while i < n:
        c = sec[i]
        if c in '0#.' or (c in 'Ee' and num is not None and i + 1 < n and sec[i + 1] in '+-'):
            if num is None:
                num = ''
                out.append(['num', ''])
            if c in 'Ee':
                out[-1][1] += sec[i:i + 2]
                i += 2
                while i < n and sec[i] in '0#':
                    out[-1][1] += sec[i]
                    i += 1
                continue
            out[-1][1] += c
        else:
            if c == '"':
                j = sec.index('"', i + 1)
                lit = sec[i + 1:j]
                i = j
            elif c == '\\' and i + 1 < n:
                lit = sec[i + 1]
                i += 1
            else:
                lit = c
            if out and out[-1][0] == 'lit':
                out[-1][1] += lit
            else:
                out.append(['lit', lit])
            if num is not None:
                num = 'done'
        i += 1
    return out


def _fmt_num(d, spec):
    """format the non-negative Decimal d with a numeric spec"""
    exp_part = None
    for k, c in enumerate(spec):
        if c in 'Ee':
            exp_part = spec[k:]
            spec = spec[:k]
            break
    if '.' in spec:
        ip, fp = spec.split('.', 1)
        has_point = True
    else:
        ip, fp = spec, ''
        has_point = False
    nint = len(ip)
    nfrac = len(fp)
    if exp_part is not None:
        esign = exp_part[1]
        edig = len(exp_part) - 2
        if d == 0:
            e = 0
            m = Decimal(0)
        else:
            e = d.adjusted() - (nint - 1 if nint > 0 else -1)
            m = _round(d.scaleb(-e), nfrac)
            lim = Decimal(10) ** (nint if nint > 0 else 0)
            if m >= lim:
                e += 1
                m = _round(d.scaleb(-e), nfrac)
        ms = _mantissa(m, '0' * len(ip), fp, has_point)
        es = '%0*d' % (edig, abs(e))
        if e < 0:
            es = '-' + es
        elif esign == '+':
            es = '+' + es
        return ms + exp_part[0] + es
    m = _round(d, nfrac)
    return _mantissa(m, ip, fp, has_point)


def _mantissa(m, ip, fp, has_point):
    nfrac = len(fp)
    s = format(m, 'f')
    if '.' in s:
        si, sf = s.split('.')
    else:
        si, sf = s, ''
    sf = (sf + '0' * nfrac)[:nfrac]
    if si == '0':
        si = ''
    # integer part
    nmin = ip.count('0')
    if len(si) < nmin:
        si = '0' * (nmin - len(si)) + si
    # fractional part: trailing '#' are optional digits
    if nfrac:
        keep = nfrac
        while keep > 0 and fp[keep - 1] == '#' and sf[keep - 1] == '0':
            keep -= 1
        sf = sf[:keep]
    res = si
    if has_point:
        res += '.'
        res += sf
    return res


def vb_format(x, fmt, single=False):
    """Format$(x, fmt) for numeric formats.  Unlike Str$/CStr, Format$
    works on the exact binary value of a Single (131.904! gives 131.90401
    with '###.#####'), so ``single`` makes no difference here."""
    secs = _split_sections(fmt)
    if isinstance(x, Decimal):
        d = x
    elif isinstance(x, int) and not isinstance(x, bool):
        d = Decimal(x)
    else:
        d = _digits15(float(x))
    neg = d < 0
    if neg and len(secs) > 1 and secs[1] != '':
        sec = secs[1]
        explicit_sign = True
    else:
        sec = secs[0]
        explicit_sign = False
    if d == 0 and len(secs) > 2 and secs[2] != '':
        sec = secs[2]
    parts = _parse(sec)
    out = ''
    ad = abs(d)
    for kind, text in parts:
        if kind == 'lit':
            out += text
        else:
            out += _fmt_num(ad, text)
    if neg and not explicit_sign:
        # VB puts the minus sign before the first digit placeholder
        out = '-' + out
    return out
