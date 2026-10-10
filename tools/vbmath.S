        .intel_syntax noprefix
        .text
# All functions: double args in xmm0 (and xmm1), result in xmm0 (SysV and Win64).
# Frame: [rsp+0] scratch double, [rsp+8] scratch double, [rsp+16] saved CW, [rsp+18] new CW
        .globl vb_pow
vb_pow:                                 # x > 0 (caller handles the rest)
        sub rsp, 40
        fnstcw word ptr [rsp+16]
        mov word ptr [rsp+18], 0x133f
        fldcw word ptr [rsp+18]
        movsd qword ptr [rsp], xmm1
        movsd qword ptr [rsp+8], xmm0
        fld qword ptr [rsp]
        fld qword ptr [rsp+8]
        fyl2x
        jmp core
        .globl vb_exp
vb_exp:
        sub rsp, 40
        fnstcw word ptr [rsp+16]
        mov word ptr [rsp+18], 0x133f
        fldcw word ptr [rsp+18]
        movsd qword ptr [rsp], xmm0
        fld qword ptr [rsp]
        fldl2e
        fmulp st(1), st
core:                                   # st0 = z ; result 2^z
        fld st(0)
        frndint                         # n, z
        fxch st(1)                      # z, n
        fsub st, st(1)                  # f = z - n, n
        ftst
        fnstsw ax
        fabs
        f2xm1
        fld1
        faddp st(1), st                 # p = 1 + (2^|f| - 1)
        test ah, 1
        jz 1f
        fld1
        fdivrp st(1), st                # p = 1 / p
1:      fscale                          # p * 2^n
        fstp st(1)                      # drop n
        fstp qword ptr [rsp]
        fldcw word ptr [rsp+16]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_log
vb_log:
        sub rsp, 40
        movsd qword ptr [rsp], xmm0
        fldln2
        fld qword ptr [rsp]
        fyl2x
        fstp qword ptr [rsp]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_logdiv                # Log(x) / c with 53-bit precision control
vb_logdiv:
        sub rsp, 40
        fnstcw word ptr [rsp+16]
        mov word ptr [rsp+18], 0x137f
        fldcw word ptr [rsp+18]
        movsd qword ptr [rsp], xmm0
        movsd qword ptr [rsp+8], xmm1
        fldln2
        fld qword ptr [rsp]
        fyl2x
        fdiv qword ptr [rsp+8]
        fstp qword ptr [rsp]
        fldcw word ptr [rsp+16]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_atn
vb_atn:
        sub rsp, 40
        movsd qword ptr [rsp], xmm0
        fld qword ptr [rsp]
        fld1
        fpatan
        fstp qword ptr [rsp]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_sin
vb_sin:
        sub rsp, 40
        movsd qword ptr [rsp], xmm0
        fld qword ptr [rsp]
        fsin
        fstp qword ptr [rsp]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_cos
vb_cos:
        sub rsp, 40
        movsd qword ptr [rsp], xmm0
        fld qword ptr [rsp]
        fcos
        fstp qword ptr [rsp]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
        .globl vb_sqrt                  # Sqr with 64-bit precision control (TRIM: 0x137f)
vb_sqrt:
        sub rsp, 40
        fnstcw word ptr [rsp+16]
        mov word ptr [rsp+18], 0x137f
        fldcw word ptr [rsp+18]
        movsd qword ptr [rsp], xmm0
        fld qword ptr [rsp]
        fsqrt
        fstp qword ptr [rsp]
        fldcw word ptr [rsp+16]
        movsd xmm0, qword ptr [rsp]
        add rsp, 40
        ret
