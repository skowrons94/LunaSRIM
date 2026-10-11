"""Math functions of the Visual Basic 5 runtime (msvbvm50.dll).

SRIM is a VB5 program.  Its ``Exp``, ``Log``, ``Atn``, ``Sin``, ``Cos`` and
``^`` do not give the correctly rounded results of the C library: they are
x87 instruction sequences (Microsoft C runtime), evaluated with 64-bit
precision and rounded to a Double at the end.  The results differ from
Python's ``math`` module in the last bit for a fraction of a percent of the
arguments, and in a Monte Carlo calculation one such bit can change a
decision (an iteration count, a flight path cut at a layer boundary) and the
rest of a trajectory.

On x86-64 the same instruction sequences are run natively (machine code
called through ctypes).  The x87 transcendental instructions are not
correctly rounded and differ in the last bit between CPU vendors; the
native routines are used only when they reproduce SRIM's values on an Intel
CPU (where the reference outputs were made).  Otherwise (or if executable
memory cannot be allocated) the functions fall back to the ``math`` module,
and the square root to an exact 64-bit-precision emulation.

The routines (msvbvm50.dll, SRIM-2013 setup):

* ``x ^ y`` (``__vbaPowerR8``, ``__vbaVarPow``) and ``Exp``: control word
  0x133f; z = y*log2(x) (``fyl2x``) or x*log2(e); n = frndint(z);
  p = 1 + f2xm1(|z - n|), inverted when z - n < 0; fscale; stored as Double.
* ``Log``: ``fldln2; fyl2x``; ``Atn``: ``fld1; fpatan``; ``Sin``/``Cos``:
  ``fsin``/``fcos``; the 80-bit result is stored to a Double by TRIM, except
  for ``Log(x) / c`` (energy grid, Log10) which is divided first.
* ``Sqr``: ``fsqrt``.  TRIM runs with the x87 control word 0x137f (64-bit
  precision), so the square root is rounded to 64 bits and then to a
  Double: for arguments just below 1 (normalised direction cosines) this
  often differs from the correctly rounded ``math.sqrt``.
"""
import math
import os
import platform

__all__ = ['pow', 'exp', 'log', 'log_div', 'atn', 'sin', 'cos', 'sqrt', 'NATIVE']


def _two_square(r):
    """r*r = hi + lo exactly (Dekker)"""
    hi = r * r
    c = 134217729.0 * r            # 2^27 + 1
    a1 = c - (c - r)
    a2 = r - a1
    lo = ((a1 * a1 - hi) + 2 * a1 * a2) + a2 * a2
    return hi, lo


def _sqrt64(x):
    """sqrt rounded to 64 bits and then to 53 bits (x87, 64-bit precision)"""
    r = math.sqrt(x)
    if not (0 < x < math.inf) or r < 1e-150 or r > 1e150:
        return r
    hi, lo = _two_square(r)
    d = (x - hi) - lo                       # x - r*r (exact enough near ties)
    if d == 0:
        return r
    ulp = math.ulp(r)
    t = abs(d) / (r + r)                    # |sqrt(x) - r|
    if abs(t - ulp / 2) > ulp / 256:
        return r
    # exact test: does sqrt(x) round (64 bits) to the midpoint r +- ulp/2 ?
    from fractions import Fraction
    sg = 1 if d > 0 else -1
    mid = Fraction(r) + sg * Fraction(ulp) / 2
    g = Fraction(ulp) / 4096
    X = Fraction(x)
    if (mid - g) ** 2 <= X <= (mid + g) ** 2:
        other = r + sg * ulp
        # ties to even (53-bit mantissa)
        return r if (math.frexp(r)[0] * 2 ** 53) % 2 == 0 else other
    return r


# x86-64 machine code (tools/vbmath.S); offsets of the entry points
_CODE = bytes.fromhex(
    '4883ec28d97c241066c74424123f13d96c2412f20f110c24f20f11442408dd0424dd442408d9f1eb'
    '1f4883ec28d97c241066c74424123f13d96c2412f20f110424dd0424d9eadec9d9c0d9fcd9c9d8e1'
    'd9e4dfe0d9e1d9f0d9e8dec1f6c4017404d9e8def1d9fdddd9dd1c24d96c2410f20f1004244883c4'
    '28c34883ec28f20f110424d9eddd0424d9f1dd1c24f20f1004244883c428c34883ec28d97c241066'
    'c74424127f13d96c2412f20f110424f20f114c2408d9eddd0424d9f1dc742408dd1c24d96c2410f2'
    '0f1004244883c428c34883ec28f20f110424dd0424d9e8d9f3dd1c24f20f1004244883c428c34883'
    'ec28f20f110424dd0424d9fedd1c24f20f1004244883c428c34883ec28f20f110424dd0424d9ffdd'
    '1c24f20f1004244883c428c34883ec28d97c241066c74424127f13d96c2412f20f110424dd0424d9'
    'fadd1c24d96c2410f20f1004244883c428c3')
_ENTRY = {'pow': 0, 'exp': 41, 'log': 122, 'log_div': 151, 'atn': 209, 'sin': 238, 'cos': 265,
          'sqrt': 292}

_keep = []          # executable buffer, kept alive

# values given by msvbvm50.dll (and by the routines below) on an Intel CPU
_INTEL_X87 = [
    ('pow', (6432.957985524626, 2.0), 41382948.44352505),
    ('pow', (145265.48020290735, 2.0), 21102059738.581264),
    ('pow', (2.6329444762163088e-05, 2.0), 6.932396614837973e-10),
    ('atn', (0.09014091311591388,), 0.08989795332127876),
]


def _load_native():
    if os.environ.get('LUNASRIM_NO_NATIVE_MATH'):
        return None
    if platform.machine().lower() not in ('x86_64', 'amd64'):
        return None
    import ctypes
    if ctypes.sizeof(ctypes.c_void_p) != 8:
        return None
    size = 4096
    if os.name == 'nt':
        k32 = ctypes.windll.kernel32
        k32.VirtualAlloc.restype = ctypes.c_void_p
        k32.VirtualAlloc.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_ulong, ctypes.c_ulong]
        addr = k32.VirtualAlloc(None, size, 0x3000, 0x40)     # MEM_COMMIT|RESERVE, RWX
        if not addr:
            return None
        ctypes.memmove(addr, _CODE, len(_CODE))
    else:
        import mmap
        buf = mmap.mmap(-1, size, prot=mmap.PROT_READ | mmap.PROT_WRITE | mmap.PROT_EXEC)
        buf.write(_CODE)
        _keep.append(buf)
        addr = ctypes.addressof(ctypes.c_char.from_buffer(buf))
    f1 = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    f2 = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    fns = {}
    for name, off in _ENTRY.items():
        proto = f2 if name in ('pow', 'log_div') else f1
        fns[name] = proto(addr + off)
    # self test: the routines must give sensible values
    if not (fns['pow'](2.0, 10.0) == 1024.0 and fns['exp'](0.0) == 1.0 and
            abs(fns['log'](10.0) - math.log(10.0)) < 1e-15 and
            abs(fns['atn'](1.0) - math.pi / 4) < 1e-15 and
            abs(fns['sin'](1.0) - math.sin(1.0)) < 1e-15 and
            abs(fns['cos'](1.0) - math.cos(1.0)) < 1e-15 and
            abs(fns['log_div'](10.0, 2.0) - math.log(10.0) / 2) < 1e-15 and
            fns['sqrt'](4.0) == 2.0 and fns['sqrt'](1 - 2 ** -53) == 1.0 and
            fns['pow'](3.0, -2.0) == math.pow(3.0, -2.0)):
        return None
    # The x87 transcendental instructions (fyl2x, f2xm1, fpatan, ...) are not
    # correctly rounded and their last bits differ between CPU vendors.  The
    # reference outputs were made on an Intel CPU; elsewhere only the square
    # root (exactly specified) is used natively.
    if not all(fns[name](*args) == v for name, args, v in _INTEL_X87):
        return {'sqrt': fns['sqrt']}
    return fns


try:
    _native = _load_native()
except Exception:       # no executable memory, ctypes missing, ...
    _native = None

NATIVE = _native is not None and 'pow' in _native
NATIVE_SQRT = _native is not None

if NATIVE:
    _npow = _native['pow']
    exp = _native['exp']
    sqrt = _native['sqrt']
    log = _native['log']
    log_div = _native['log_div']
    atn = _native['atn']
    sin = _native['sin']
    cos = _native['cos']

    def pow(x, y):
        """VB ``x ^ y``"""
        x = float(x)
        y = float(y)
        if x > 0:
            if y == 0:
                return 1.0
            if math.isinf(x) or math.isinf(y) or x != x or y != y:
                return math.pow(x, y)
            return _npow(x, y)
        if x == 0 or y == 0 or y != math.floor(y):
            return math.pow(x, y)
        # negative base, integer exponent: |x| ^ y with the sign of the result
        r = _npow(-x, y)
        return -r if math.fmod(y, 2) else r
else:
    sqrt = _native['sqrt'] if NATIVE_SQRT else _sqrt64
    exp = math.exp
    log = math.log
    atn = math.atan
    sin = math.sin
    cos = math.cos
    pow = math.pow

    def log_div(x, c):
        return math.log(x) / c
