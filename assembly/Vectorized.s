#define STDOUT 0xd0580000

.section .data

real:           .float 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0
imag:           .float 0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 1.0
twiddle_real:   .float 1.0,  0.7071,  0.0, -0.7071,  -1.0, -0.7071,  0.0,  0.7071
twiddle_imag:   .float 0.0, -0.7071, -1.0, -0.7071,   0.0,  0.7071,  1.0,  0.7071
bitrev_indices: .space 32          # 8 words (for n = 8)

filename_vl_debug:   .string "vl_output.hex"
.byte 0,0               # pad to 16 bytes

.align 4
vl_debug_value:   .word 0xA5A5A5A5

.align 4
filename_rev_real0: .string "rev_real0_output.hex"
.byte 0,0,0            # pad to 24 bytes

.align 4
filename_real:    .string "output_real.hex"
.byte 0,0,0           # pad to 20 bytes

.align 4
filename_imag:    .string "output_imag.hex"
.byte 0,0,0           # pad to 20 bytes

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

fft_stage_loop:
    li   t0, 8                # N = 8
    li   t6, 3                # log2(N) = 3
    bgt  s0, t6, next_stage   # Skip if stage > log2(N)
    
    sll  s1, s0, 1            # m = 2^stage
    bgt  s1, t0, next_stage   # Skip if m > N
    
    srli s2, s1, 1            # half_m = m/2
    div  s3, t0, s1           # stride = N/m
    li   s4, 0                # i = 0

fft_i_loop:
    bge  s4, t0, next_stage   
    li   s5, 0                # j = 0

fft_j_loop:
    bge  s5, s2, fft_i_inc    
    
    # Calculate indices with bounds checking
    add  t1, s4, s5           # i+j
    bge  t1, s11, fft_j_inc   
    
    add  t2, t1, s2           # i+j+half_m
    bge  t2, s11, fft_j_inc   
    
    slli t3, t1, 2            # byte offset top
    slli t4, t2, 2            # byte offset bottom

    # Load elements using vector loads
    add  t5, a0, t3
    vle32.v v0, (t5)          # real[i+j]
    add  t5, a1, t3
    vle32.v v1, (t5)          # imag[i+j]

    add  t5, a0, t4
    vle32.v v2, (t5)          # real[i+j+half_m]
    add  t5, a1, t4
    vle32.v v3, (t5)          # imag[i+j+half_m]

    # Load twiddle with bounds check
    mul  t5, s5, s3           
    bge  t5, s11, fft_j_inc   
    slli t5, t5, 2            
    add  t6, a2, t5
    vle32.v v4, (t6)          # W_real
    add  t6, a3, t5
    vle32.v v5, (t6)          # W_imag

    # Store to temp scalar registers
    vfmv.f.s ft0, v0          # real_t
    vfmv.f.s ft1, v1          # imag_t
    vfmv.f.s ft2, v2          # real_b
    vfmv.f.s ft3, v3          # imag_b
    vfmv.f.s ft4, v4          # W_real
    vfmv.f.s ft5, v5          # W_imag

    # Butterfly computation in scalar
    fmul.s ft6, ft2, ft4      # real_b * W_real
    fmul.s ft7, ft3, ft5      # imag_b * W_imag
    fsub.s ft8, ft6, ft7      # temp_real
    fmul.s ft6, ft2, ft5      # real_b * W_imag
    fmul.s ft7, ft3, ft4      # imag_b * W_real
    fadd.s ft9, ft6, ft7      # temp_imag

    # Move results back to vector registers
    fsub.s ft10, ft0, ft8     # real_b_new
    vfmv.s.f v12, ft10
    add  t5, a0, t4
    vse32.v v12, (t5)

    fsub.s ft11, ft1, ft9     # imag_b_new
    vfmv.s.f v13, ft11
    add  t5, a1, t4
    vse32.v v13, (t5)

    fadd.s ft10, ft0, ft8     # real_t_new
    vfmv.s.f v10, ft10
    add  t5, a0, t3
    vse32.v v10, (t5)

    fadd.s ft11, ft1, ft9     # imag_t_new
    vfmv.s.f v11, ft11
    add  t5, a1, t3
    vse32.v v11, (t5)

fft_j_inc:
    addi s5, s5, 1            
    j    fft_j_loop

fft_i_inc:
    add  s4, s4, s1           
    j    fft_i_loop

next_stage:
    addi s0, s0, 1            
    ble  s0, t6, fft_stage_loop

    # Step 5: Write output files
    la   a0, filename_real
    la   a1, real
    li   a2, 32               # 8 floats * 4 bytes
    call write_to_file

    la   a0, filename_imag
    la   a1, imag
    li   a2, 32
    call write_to_file

    j    _finish

write_to_file:
    # Instead of file operations, write directly to memory
    addi sp, sp, -16
    sw   ra, 12(sp)
    sw   s0, 8(sp)
    sw   s1, 4(sp)
    sw   s2, 0(sp)

    mv   s0, a0    # Save output address
    mv   s1, a1    # Save data address
    mv   s2, a2    # Save size

    # Direct memory copy loop
    mv   t0, zero  # Initialize counter
copy_loop:
    bge  t0, s2, copy_done
    add  t1, s1, t0          # Source address
    flw  ft0, 0(t1)          # Load float
    add  t2, s0, t0          # Destination address
    fsw  ft0, 0(t2)          # Store float
    addi t0, t0, 4           # Increment by float size
    j    copy_loop

copy_done:
    lw   ra, 12(sp)
    lw   s0, 8(sp)
    lw   s1, 4(sp)
    lw   s2, 0(sp)
    addi sp, sp, 16
    ret

_finish:
    lui  t0, 0xd0580          # STDOUT address
    li   t1, 3                # Success code
    sb   t1, 0(t0)            # Write to STDOUT

_final_halt_loop:
    j    _final_halt_loop