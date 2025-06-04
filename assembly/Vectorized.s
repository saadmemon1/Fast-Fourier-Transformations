#define STDOUT 0xd0580000

.section .data

real:           .float 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0
imag:           .float 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
twiddle_real:   .float 1.0, 1.0, 0.0, 1.0, 0.7071, 0.0, -0.7071
twiddle_imag:   .float 0.0, 0.0, -1.0, 0.0, -0.7071, -1.0, -0.7071
bitrev_indices: .space 32          # 8 words (for n = 8)

.align 4
realhex:    .string "./veer/tempFiles/real.hex"
.byte 0,0,0

.align 4
imaghex:    .string "./veer/tempFiles/imag.hex"
.byte 0,0,0

.section .text.init
.globl _start

_start:
    li   s11, 8              # n = 8

    # Step 1: Compute bit-reversal indices
    la   a2, bitrev_indices
    li   t0, 0
bitrev_loop:
    mv   t3, t0
    li   t2, 0
    li   t6, 3             # number of bits = 3
reverse_bits:
    slli t2, t2, 1
    andi t5, t3, 1
    or   t2, t2, t5
    srli t3, t3, 1
    addi t6, t6, -1
    bnez t6, reverse_bits
    sw   t2, 0(a2)
    addi t0, t0, 1
    addi a2, a2, 4
    blt  t0, s11, bitrev_loop

    # Step 2: Bit-reversal permutation (in-place swap)
    la   t0, real
    la   t1, imag
    la   t2, bitrev_indices
    li   t3, 0                # i = 0
bitrev_perm_loop:
    bge  t3, s11, end_bitrev  # i >= n? done
    slli t4, t3, 2            # byte offset for i
    add  t5, t2, t4
    lw   t5, 0(t5)            # j = bitrev_indices[i]
    bge  t3, t5, next_i       # skip if i >= j
    slli t6, t5, 2            # byte offset for j
    # Swap real[i] and real[j]
    add  a0, t0, t4
    add  a1, t0, t6
    flw  ft0, 0(a0)
    flw  ft1, 0(a1)
    fsw  ft0, 0(a1)
    fsw  ft1, 0(a0)
    # Swap imag[i] and imag[j]
    add  a0, t1, t4
    add  a1, t1, t6
    flw  ft0, 0(a0)
    flw  ft1, 0(a1)
    fsw  ft0, 0(a1)
    fsw  ft1, 0(a0)
next_i:
    addi t3, t3, 1
    j    bitrev_perm_loop
end_bitrev:
    # Step 3: Vector configuration
    li   t0, 1                # Vector length = 1
    vsetvli zero, t0, e32, m1, ta, ma

    # Step 4: FFT stages with bounds checking
    li   s0, 1                # stage = 1
    la   a0, real
    la   a1, imag
    la   a2, twiddle_real
    la   a3, twiddle_imag
    
    # Stage identifiers
    li   s8, 1                # Stage 1
    li   s9, 2                # Stage 2
    li   s10, 4               # Stage 3

fft_stage_loop:
    li   t0, 8                # N = 8
    li   t6, 3                # log2(N) = 3
    bgt  s0, t6, next_stage
    
    sll  s1, s0, 1            # m = 2^stage
    srli s2, s1, 1            # half_m = m/2
    li   s4, 0                # group offset = 0
    
    # Set twiddle start index
    beq  s0, s8, stage1_twiddle
    beq  s0, s9, stage2_twiddle
    beq  s0, s10, stage3_twiddle
    j    group_setup

stage1_twiddle:
    li   s7, 0                # Use W₈⁰
    j    group_setup

stage2_twiddle:
    li   s7, 0                # Use W₈⁰, W₈²
    j    group_setup

stage3_twiddle:
    li   s7, 0                # Use W₈⁰, W₈¹, W₈², W₈³

group_setup:
group_loop:
    li   s5, 0                # j = 0

butterfly_loop:
    bge  s5, s2, next_group
    
    # Calculate indices
    add  t1, s4, s5           # even_idx = group + j
    add  t2, t1, s2           # odd_idx = even_idx + half_m
    
    # Load elements
    slli t3, t1, 2            # even offset
    slli t4, t2, 2            # odd offset
    
    add  t5, a0, t3
    flw  ft0, 0(t5)           # even_real
    add  t5, a1, t3
    flw  ft1, 0(t5)           # even_imag
    
    add  t5, a0, t4
    flw  ft2, 0(t5)           # odd_real
    add  t5, a1, t4
    flw  ft3, 0(t5)           # odd_imag

    # Calculate twiddle index based on stage and position
    mul  t4, s5, t0           # j * N
    div  t4, t4, s1           # (j * N) / m
    add  t4, t4, s7           # Add stage offset
    slli t4, t4, 2            # Convert to byte offset

    # Load twiddle factors
    add  t5, a2, t4
    flw  ft4, 0(t5)           # W_real
    add  t5, a3, t4
    flw  ft5, 0(t5)           # W_imag

    # Complex multiplication
    fmul.s ft6, ft2, ft4      # odd_real * W_real
    fmul.s ft7, ft3, ft5      # odd_imag * W_imag
    fmul.s ft8, ft2, ft5      # odd_real * W_imag
    fmul.s ft9, ft3, ft4      # odd_imag * W_real
    
    fsub.s ft6, ft6, ft7      # temp_real = re*wr - im*wi
    fadd.s ft7, ft8, ft9      # temp_imag = re*wi + im*wr

    # Butterfly computation
    fadd.s ft8, ft0, ft6      # even' = even + temp
    fadd.s ft9, ft1, ft7
    fsub.s ft10, ft0, ft6     # odd' = even - temp
    fsub.s ft11, ft1, ft7

    # Store results
    slli t3, t1, 2
    add  t5, a0, t3
    fsw  ft8, 0(t5)           # store even_real
    add  t5, a1, t3
    fsw  ft9, 0(t5)           # store even_imag
    
    slli t3, t2, 2
    add  t5, a0, t3
    fsw  ft10, 0(t5)          # store odd_real
    add  t5, a1, t3
    fsw  ft11, 0(t5)          # store odd_imag

    addi s5, s5, 1            # j++
    j    butterfly_loop

next_group:
    add  s4, s4, s1           # group += m
    blt  s4, t0, group_loop

next_stage:
    slli s0, s0, 1            # stage *= 2
    ble  s0, t6, fft_stage_loop

    # Write exactly 64 bytes (8 complex values)
    li   t0, 0xd0580100      # Output memory region
    la   t1, real            # Source real
    la   t2, imag            # Source imag
    li   t3, 0               # Counter
    li   t4, 32              # 8 floats * 4 bytes

write_loop:
    bge  t3, t4, write_done
    
    # Write interleaved real/imag
    add  t5, t1, t3          # real[i] address
    add  t6, t2, t3          # imag[i] address
    flw  ft0, 0(t5)          # Load real
    flw  ft1, 0(t6)          # Load imag
    
    add  t5, t0, t3
    slli t6, t3, 1          # Double offset for complex pairs
    add  t5, t5, t6         # Output address
    
    fsw  ft0, 0(t5)         # Store real
    fsw  ft1, 4(t5)         # Store imag next to real
    
    addi t3, t3, 4
    j    write_loop

write_done:
    lui  t0, 0xd0580
    li   t1, 1
    sw   t1, 0(t0)

_final_halt_loop:
    j    _final_halt_loop

.end
