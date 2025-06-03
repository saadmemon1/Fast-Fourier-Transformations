# --- Data Section ---
.data 0x10000000
.align 2
.equ dataSize, 8
.equ halfDataSize, 4

real: 
    .float 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0   # enough room for N=8; 
                                                  # if N=4, only the first 4 floats matter
imag: 
    .float 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

bitrev_indices: 
    .space 32    # enough for up to N=8 (8×4=32 bytes)

# Mathematical constants (used nowhere except as reference)
PI:         .float 3.14159265358979323846
TWO_PI:     .float 6.28318530717958647692
NEG_TWO_PI: .float -6.28318530717958647692
HALF_PI:    .float 1.57079632679489661923

# Polynomial coefficients for cos(x) ≈ c0·x^8 + c1·x^6 + c2·x^4 + c3·x^2 + c4
cos_coeff_0: .float 2.44677067e-5
cos_coeff_1: .float -1.38877297e-3
cos_coeff_2: .float 4.16666567e-2
cos_coeff_3: .float -5.00000000e-1
cos_coeff_4: .float 1.00000000e+0

# Polynomial coefficients for sin(x) ≈ s0·x^7 + s1·x^5 + s2·x^3 + s3·x
sin_coeff_0: .float 2.86567956e-6
sin_coeff_1: .float -1.98559923e-4
sin_coeff_2: .float 8.33338592e-3
sin_coeff_3: .float -1.66666672e-1

# --- Text Section ---
.text
.globl main

main:
    # ------------------------------------------------
    # 1) Load polynomial coefficients into saved FP regs
    # ------------------------------------------------
    la   a3, cos_coeff_0
    flw  fs0, 0(a3)    # c0
    flw  fs1, 4(a3)    # c1
    flw  fs2, 8(a3)    # c2
    flw  fs3, 12(a3)   # c3
    flw  fs4, 16(a3)   # c4

    la   a3, sin_coeff_0
    flw  fs5, 0(a3)    # s0
    flw  fs6, 4(a3)    # s1
    flw  fs7, 8(a3)    # s2
    flw  fs8, 12(a3)   # s3

    # ------------------------------------------------
    # 2) Initialize N (power of two) in s11
    #    For example, set N = 8 (8 complex samples)
    # ------------------------------------------------
    li   s11, 8          # Change to 4 if you want N=4, or 16 for N=16, etc.

    # ------------------------------------------------
    # Step 1: Compute bit‐reversal indices (generic)
    # ------------------------------------------------
    la   a2, bitrev_indices
    li   t0, 0
    mv   t1, s11
    li   t4, 0
    mv   t3, s11
    addi t3, t3, -1
    beqz t3, bitrev_loop   # if s11-1 == 0 ⇒ t3=0 ⇒ goto bitrev_loop
count_bits:
    addi t4, t4, 1
    srli t3, t3, 1
    bnez t3, count_bits

bitrev_loop:
    mv   t3, t0
    li   t2, 0
    beqz t4, store_index_label
    mv   t6, t4
reverse_bits:
    slli  t2, t2, 1
    andi  t5, t3, 1
    or    t2, t2, t5
    srli  t3, t3, 1
    addi  t6, t6, -1
    bnez  t6, reverse_bits

store_index_label:
    sw    t2, 0(a2)     # store bit‐reversed index
    addi  t0, t0, 1
    addi  a2, a2, 4
    blt   t0, t1, bitrev_loop

    # ------------------------------------------------
    # Step 2: Bit‐reversal permutation (generic)
    # ------------------------------------------------
    la   a0, real       # a0 → base of real[]
    la   a1, imag       # a1 → base of imag[]
    la   a2, bitrev_indices
    li   t0, 0

permute_loop:
    lw    t1, 0(a2)     # t1 = bitrev_indices[t0]
    bge   t0, t1, skip_swap_label
    # swap real[t0] <-> real[t1]
    slli  t2, t0, 2
    add   t3, a0, t2
    slli  t4, t1, 2
    add   t5, a0, t4
    flw   ft0, 0(t3)
    flw   ft1, 0(t5)
    fsw   ft1, 0(t3)
    fsw   ft0, 0(t5)
    # swap imag[t0] <-> imag[t1]
    add   t3, a1, t2
    add   t5, a1, t4
    flw   ft0, 0(t3)
    flw   ft1, 0(t5)
    fsw   ft1, 0(t3)
    fsw   ft0, 0(t5)

skip_swap_label:
    addi  a2, a2, 4
    addi  t0, t0, 1
    blt   t0, s11, permute_loop

    # ------------------------------------------------
    # Step 3: FFT stages (generic polynomial‐twiddle)
    # ------------------------------------------------
    li    s0, 1         # s0 = size of sub‐FFT = 1 (will double each stage)
    mv    s1, s11       # s1 = N

stage_loop:
    slli  s3, s0, 1        # s3 = m = 2*s0 (current FFT “butterfly width”)
    srli  s6, s3, 1        # s6 = half = m/2
    li    s2, 0           # s2 = group_start = 0 at beginning of each stage

group_loop:
    li    t0, 0           # t0 = j (butterfly index within group)

butterfly_loop:
    # Calculate indices for “even” (t1) and “odd” (t2) elements:
    add   t1, s2, t0       # t1 = group_start + j
    add   t2, t1, s6       # t2 = t1 + half

    # Load “even” sample (real+t1, imag+t1):
    slli  t3, t1, 2
    add   t4, a0, t3
    add   t5, a1, t3
    flw   ft0, 0(t4)       # ft0 = real[t1]
    flw   ft1, 0(t5)       # ft1 = imag[t1]

    # Load “odd” sample (real+t2, imag+t2):
    slli  t3, t2, 2
    add   t6, a0, t3
    add   t3, a1, t3
    flw   ft2, 0(t6)       # ft2 = real[t2]
    flw   ft3, 0(t3)       # ft3 = imag[t2]

    # ----------------------------
    # Compute angle = –2·π·j / m
    #   Let:
    #     ft4 = (float) j
    #     ft5 = (float) m
    #     ft6 = ft4 / ft5 = j/m
    #     ft7 = TWO_PI constant in ft7 (bitpattern 0x40C90FDB)
    #     fmul: ft6 = (j/m) × (2π) → ft6 = angle
    #     fneg: ft6 = –angle
    # ----------------------------
    fcvt.s.w   ft4,  t0       # ft4 = float(j)
    fcvt.s.w   ft5,  s3       # ft5 = float(m)
    fdiv.s     ft6,  ft4, ft5 # ft6 = j/m
    li         t4,   0x40C90FDB  # 0x40C90FDB is IEEE-754 bitpattern for 2·π
    fmv.s.x    ft7,  t4         # ft7 = 2π
    fmul.s     ft6,  ft6, ft7   # ft6 = (j/m) * (2π)
    fneg.s     ft6,  ft6         # ft6 = –(2π·j/m) = angle

    # ----------------------------
    # Evaluate cosine polynomial: ft8 = cos(angle)
    #   Using coefficients in fs0..fs4
    #   Let x = ft6; x² = ft4 (reusing ft4)
    # ----------------------------
    fmul.s   ft4,  ft6, ft6    # ft4 = x²
    fmv.s    ft8,  fs0         # ft8 = c0
    fmul.s   ft8,  ft8, ft4    # c0 * x²
    fadd.s   ft8,  ft8, fs1    # + c1
    fmul.s   ft8,  ft8, ft4    # * x²
    fadd.s   ft8,  ft8, fs2    # + c2
    fmul.s   ft8,  ft8, ft4    # * x²
    fadd.s   ft8,  ft8, fs3    # + c3
    fmul.s   ft8,  ft8, ft4    # * x²
    fadd.s   ft8,  ft8, fs4    # + c4
    # → ft8 now holds W_real = cos(angle)

    # ----------------------------
    # Evaluate sine polynomial: ft5 = sin(angle)
    #   Using coefficients in fs5..fs8
    #   Reuse ft4 = x² from above
    # ----------------------------
    fmv.s    ft5,  fs5         # ft5 = s0
    fmul.s   ft5,  ft5, ft4    # s0 * x²
    fadd.s   ft5,  ft5, fs6    # + s1
    fmul.s   ft5,  ft5, ft4    # * x²
    fadd.s   ft5,  ft5, fs7    # + s2
    fmul.s   ft5,  ft5, ft4    # * x²
    fadd.s   ft5,  ft5, fs8    # + s3
    fmul.s   ft5,  ft5, ft6    # * x
    # → ft5 now holds W_imag = sin(angle)

    # ----------------------------
    # Complex multiply: (odd_real + j·odd_imag) × (W_real + j·W_imag)
    #   real_part = odd_real*W_real – odd_imag*W_imag
    #   imag_part = odd_real*W_imag + odd_imag*W_real
    # ----------------------------
    fmul.s   ft6,  ft2,  ft8    # ft6 = odd_real * W_real
    fmul.s   ft7,  ft3,  ft5    # ft7 = odd_imag * W_imag
    fsub.s   ft6,  ft6,  ft7    # ft6 = real_part

    fmul.s   ft7,  ft2,  ft5    # ft7 = odd_real * W_imag
    fmul.s   ft9,  ft3,  ft8    # ft9 = odd_imag * W_real
    fadd.s   ft7,  ft7,  ft9    # ft7 = imag_part

    # ----------------------------
    # Butterfly: combine “even” (ft0,ft1) plus/minus (ft6,ft7):
    #   new_even_real = ft0 + ft6
    #   new_even_imag = ft1 + ft7
    #   new_odd_real  = ft0 – ft6
    #   new_odd_imag  = ft1 – ft7
    # ----------------------------
    fadd.s  ft9,   ft0, ft6     # ft9  = new_even_real
    fadd.s  ft10,  ft1, ft7     # ft10 = new_even_imag
    fsub.s  ft11,  ft0, ft6     # ft11 = new_odd_real
    fsub.s  ft3,   ft1, ft7     # ft3  = new_odd_imag

    # ----------------------------
    # Store results back to memory:
    #   real[t1] ← ft9, imag[t1] ← ft10
    #   real[t2] ← ft11, imag[t2] ← ft3
    # ----------------------------
    slli   t3,  t1,  2
    add    t4,  a0,  t3
    add    t5,  a1,  t3
    fsw    ft9,  0(t4)
    fsw    ft10, 0(t5)

    slli   t3,  t2,  2
    add    t4,  a0,  t3
    add    t5,  a1,  t3
    fsw    ft11, 0(t4)
    fsw    ft3,  0(t5)

    # ----------------------------
    # Next butterfly index j += 1
    # ----------------------------
    addi   t0,  t0,  1
    blt    t0,  s6,  butterfly_loop

    # ----------------------------
    # Next group_start += m
    # ----------------------------
    add    s2,  s2,  s3
    blt    s2,  s1,  group_loop

    # ----------------------------
    # Next stage: m doubles (s0 <<= 1)
    # ----------------------------
    slli   s0,  s0,  1
    blt    s0,  s1,  stage_loop

done:
    j done