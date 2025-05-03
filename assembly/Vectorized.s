.section .text
.global fft_stage_vectorized
fft_stage_vectorized:
    # a0 = base address of matrix (interleaved complex)
    # a1 = stage size (s0)
    # a2 = FFT length N
    # a3 = address of twiddle_re
    # a4 = address of twiddle_im

    # Constants
    slli t0, a1, 1              # t0 = stage_size = 2 * s0
    li t1, 0                    # block base index

stage_block_loop:
    bge t1, a2, stage_done      # if block_start >= N, done

    li t2, 0                    # within-block offset k
stage_inner_loop:
    bge t2, a1, next_block      # if k >= s0, move to next block

    # index = block_base + k
    add t3, t1, t2
    slli t4, t3, 3              # offset in bytes (8 bytes per complex)
    add t5, a0, t4              # addr1 = &x[k]

    # index2 = index + s0
    add t6, t3, a1
    slli t7, t6, 3
    add t8, a0, t7              # addr2 = &x[k + s0]

    # Load x[k] and x[k+s0] as vectors (Re, Im pairs)
    vsetvli t9, zero, e32, m1
    vle32.v v0, (t5)            # v0 = [Re_even, Im_even]
    vle32.v v1, (t8)            # v1 = [Re_odd, Im_odd]

    # Load twiddle factor for this k
    mul t10, t2, 4              # offset = 4 * k
    add t11, a3, t10
    flw f0, 0(t11)              # f0 = cos
    add t11, a4, t10
    flw f1, 0(t11)              # f1 = -sin

    # Compute: temp = v1 * twiddle (complex multiply)
    # temp_re = Re*v0 - Im*v1
    # temp_im = Re*v1 + Im*v0

    vmv.v.x v2, f0              # v2 = [cos, cos]
    vmv.v.x v3, f1              # v3 = [-sin, -sin]

    vmul.vv v4, v1, v2          # v4 = v1 * cos
    vmul.vv v5, v1, v3          # v5 = v1 * -sin

    vslideup.vi v6, v1, 1       # shift to align real/imag for complex ops
    vmul.vv v6, v6, v3          # v6 = Re*sin or Im*sin
    vmul.vv v7, v6, v2          # v7 = Re*cos or Im*cos

    # Do fused multiply/add/subtract if supported
    # Add/sub with original even vector
    vadd.vv v8, v0, v4          # out1 = even + temp
    vsub.vv v9, v0, v4          # out2 = even - temp

    vse32.v v8, (t5)            # store x[k]
    vse32.v v9, (t8)            # store x[k+s0]

    addi t2, t2, 1
    j stage_inner_loop

next_block:
    add t1, t1, t0              # block += 2*s0
    j stage_block_loop

stage_done:
    ret
