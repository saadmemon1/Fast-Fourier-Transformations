#line no: 188, 111, 104, 54, 49, 19, 17, 14, 11, 9, 6, 5 change all these values to change n


# --- Data Section ---
.data
.align 2
.equ dataSize, 1024  # ← Change this to desired N (must be power of 2)
.equ halfDataSize, 512  # ← Change to N/2

real: 
  .space 4096             # ← Change to N*4 (e.g., 256→1024, 512→2048)
imag: 
  .space 4096               # ← Same as real

bitrev_indices: 
  .space 4096               # ← Same as real

W_real:
  .space 2048              # ← Change to (N/2)*4
W_imag:
  .space 2048              # ← Change to (N/2)*4

# Mathematical constants
PI: .float 3.14159265358979323846
TWO_PI: .float 6.28318530717958647692
NEG_TWO_PI: .float -6.28318530717958647692

# Polynomial coefficients
cos_coeff:
  .float 2.44677067e-5    # c0
  .float -1.38877297e-3   # c1
  .float 4.16666567e-2    # c2
  .float -5.00000000e-1   # c3
  .float 1.00000000e+0    # c4

sin_coeff:
  .float 2.86567956e-6    # s0
  .float -1.98559923e-4   # s1
  .float 8.33338592e-3    # s2
  .float -1.66666672e-1   # s3

# --- Text Section ---
.text
.globl main

main:
  # ================================================
  # Step 1: Compute bit-reversal indices
  # ================================================
  la a2, bitrev_indices
  li t0, 0                 # i = 0
  li t1, 1024              # N = 1024 # ← Change to N

bitrev_compute_loop:
  mv t3, t0                # Copy i
  li t2, 0                 # reversed_i = 0
  li t4, 10                # 10 bits for 1024 # ← Change to log2(N) (e.g., 8 for 256)
  
reverse_bits_loop:
  slli t2, t2, 1
  andi t5, t3, 1
  or t2, t2, t5
  srli t3, t3, 1
  addi t4, t4, -1
  bnez t4, reverse_bits_loop
  
  sw t2, 0(a2)
  addi t0, t0, 1
  addi a2, a2, 4
  blt t0, t1, bitrev_compute_loop

  # ================================================
  # Step 2: Bit-reversal permutation
  # ================================================
  la a0, real
  la a1, imag
  la a2, bitrev_indices
  li t0, 0                 # i = 0

bitrev_permute_loop:
  lw t1, 0(a2)             # reversed index
  bge t0, t1, skip_swap

  # Calculate addresses using integer registers
  slli t2, t0, 2          # i * 4
  slli t3, t1, 2          # reversed * 4
  
  # Swap real elements
  add t4, a0, t2           # real[i]
  add t5, a0, t3           # real[reversed]
  flw ft0, 0(t4)
  flw ft1, 0(t5)
  fsw ft1, 0(t4)
  fsw ft0, 0(t5)

  # Swap imag elements
  add t4, a1, t2           # imag[i]
  add t5, a1, t3           # imag[reversed]
  flw ft0, 0(t4)
  flw ft1, 0(t5)
  fsw ft1, 0(t4)
  fsw ft0, 0(t5)

skip_swap:
  addi a2, a2, 4
  addi t0, t0, 1
  li t6, 1024   # ← Change to N
  blt t0, t6, bitrev_permute_loop

  # ================================================
  # Step 3: FFT stages with dynamic twiddle calculation
  # ================================================
  li s0, 1                 # stage = 1
  li s1, 1024              # N = 1024  # ← Change to N

stage_loop:
  slli s3, s0, 1           # m = 2*stage
  srli s6, s3, 1           # m/2
  
  # Precompute twiddle factors
  mv a0, s3
  jal ra, precompute_twiddles

  li s2, 0                 # group offset

group_loop:
  li t0, 0                 # j = 0

butterfly_loop:
  add t1, s2, t0           # even index
  add t2, t1, s6           # odd index

  # Load even elements
  slli t3, t1, 2
  add t4, a0, t3           # real[even]
  add t5, a1, t3           # imag[even]
  flw ft0, 0(t4)
  flw ft1, 0(t5)

  # Load odd elements
  slli t3, t2, 2
  add t4, a0, t3           # real[odd]
  add t5, a1, t3           # imag[odd]
  flw ft2, 0(t4)
  flw ft3, 0(t5)

  # Load twiddle factors
  slli t3, t0, 2
  la a3, W_real
  add a3, a3, t3
  flw ft4, 0(a3)           # W_real
  la a3, W_imag
  add a3, a3, t3
  flw ft5, 0(a3)           # W_imag

  # Complex multiplication
  fmul.s ft6, ft2, ft4     # real part
  fmul.s ft7, ft3, ft5
  fsub.s ft6, ft6, ft7

  fmul.s ft7, ft2, ft5     # imag part
  fmul.s ft8, ft3, ft4
  fadd.s ft7, ft7, ft8

  # Update elements
  fadd.s ft8, ft0, ft6     # new_even_real
  fadd.s ft9, ft1, ft7     # new_even_imag
  fsub.s ft10, ft0, ft6    # new_odd_real
  fsub.s ft11, ft1, ft7    # new_odd_imag

  # Store results
  slli t3, t1, 2
  add t4, a0, t3
  add t5, a1, t3
  fsw ft8, 0(t4)
  fsw ft9, 0(t5)

  slli t3, t2, 2
  add t4, a0, t3
  add t5, a1, t3
  fsw ft10, 0(t4)
  fsw ft11, 0(t5)

  addi t0, t0, 1
  blt t0, s6, butterfly_loop

  add s2, s2, s3
  blt s2, s1, group_loop

  slli s0, s0, 1
  li t6, 1024   # ← Change to N
  blt s0, t6, stage_loop

  j done

# --- Twiddle Factor Precomputation ---
precompute_twiddles:
  la t0, W_real
  la t1, W_imag
  li t2, 0                 # k counter
  
  # Load constants using integer base registers
  la a4, NEG_TWO_PI
  flw ft0, 0(a4)           # -2π
  fcvt.s.w ft1, a0         # m as float
  fdiv.s ft0, ft0, ft1     # -2π/m

twiddle_loop:
  fcvt.s.w ft2, t2         # k as float
  fmul.s ft3, ft0, ft2     # angle = -2πk/m
  
  # Calculate cos(angle)
  fmv.s fa0, ft3
  jal ra, cos_approx
  fsw fa0, 0(t0)
  
  # Calculate sin(angle)
  fmv.s fa0, ft3
  jal ra, sin_approx
  fsw fa0, 0(t1)
  
  addi t2, t2, 1
  addi t0, t0, 4
  addi t1, t1, 4
  srli t3, a0, 1           # m/2
  blt t2, t3, twiddle_loop
  ret

# --- Cosine Approximation ---
cos_approx:
  # Range reduction to [-π, π]
  la a4, TWO_PI
  flw ft0, 0(a4)
  la a4, PI
  flw ft1, 0(a4)
  
  fabs.s ft2, fa0
  flt.s t0, ft2, ft1
  bnez t0, cos_poly
  
  # x = x - floor((x + π)/(2π)) * 2π
  la a4, PI
  flw ft3, 0(a4)
  fadd.s ft4, fa0, ft3
  fdiv.s ft4, ft4, ft0
  fcvt.w.s t1, ft4
  fcvt.s.w ft4, t1
  fmul.s ft4, ft4, ft0
  fsub.s fa0, fa0, ft4

cos_poly:
  # Load coefficients from array
  la a4, cos_coeff
  flw ft6, 0(a4)           # c0
  flw ft7, 4(a4)           # c1
  flw ft8, 8(a4)           # c2
  flw ft9, 12(a4)          # c3
  flw ft10, 16(a4)         # c4
  
  fmul.s ft5, fa0, fa0     # x²
  fmul.s ft6, ft6, ft5     # c0*x²
  fadd.s ft6, ft6, ft7     # c1 + c0*x²
  fmul.s ft6, ft6, ft5     # *x²
  fadd.s ft6, ft6, ft8     # c2 + ...
  fmul.s ft6, ft6, ft5     # *x²
  fadd.s ft6, ft6, ft9     # c3 + ...
  fmul.s ft6, ft6, ft5     # *x²
  fadd.s fa0, ft6, ft10    # c4 + ...
  ret

# --- Sine Approximation ---
sin_approx:
  # Range reduction similar to cosine
  la a4, TWO_PI
  flw ft0, 0(a4)
  la a4, PI
  flw ft1, 0(a4)
  
  fabs.s ft2, fa0
  flt.s t0, ft2, ft1
  bnez t0, sin_poly
  
  la a4, PI
  flw ft3, 0(a4)
  fadd.s ft4, fa0, ft3
  fdiv.s ft4, ft4, ft0
  fcvt.w.s t1, ft4
  fcvt.s.w ft4, t1
  fmul.s ft4, ft4, ft0
  fsub.s fa0, fa0, ft4

sin_poly:
  # Load coefficients from array
  la a4, sin_coeff
  flw ft6, 0(a4)           # s0
  flw ft7, 4(a4)           # s1
  flw ft8, 8(a4)           # s2
  flw ft9, 12(a4)          # s3
  
  fmul.s ft5, fa0, fa0     # x²
  fmul.s ft6, ft6, ft5     # s0*x²
  fadd.s ft6, ft6, ft7     # s1 + s0*x²
  fmul.s ft6, ft6, ft5     # *x²
  fadd.s ft6, ft6, ft8     # s2 + ...
  fmul.s ft6, ft6, ft5     # *x²
  fadd.s ft6, ft6, ft9     # s3 + ...
  fmul.s fa0, ft6, fa0     # *x
  ret

done:
  j done
