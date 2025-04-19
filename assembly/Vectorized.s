#define STDOUT 0xd0580000

.section .text
.global _start
_start:
## START YOUR CODE HERE

la a0, matrix
lw a1, size
call printToLogVectorized

la a0, matrix
lw a1, size
call transpose

la a0, matrix
lw a1, size
call printToLogVectorized


j _finish

transpose:
    # a0 = base address, a1 = n (matrix dimension)

    # Initialize outer loop index: i = 0
    addi    s0, zero, 0         # s0 holds i

outer_loop:
    bge     s0, a1, done         # if i >= n, we are done
    addi    s1, s0, 1           # initialize inner loop index: j = i + 1

inner_loop:
    bge     s1, a1, next_i       # if j >= n, advance to next row

    # Compute remaining elements count: t1 = n - j
    sub     t1, a1, s1

    # Set vector length (VL) for 32-bit floats; t0 = actual vl = min(n-j, VLEN)
    vsetvli t0, t1, e32, m1

    # --- Compute address for matrix[i][j] -----------------------
    # Element (i,j) address = a0 + ((i*n + j) * 4)
    mul     t2, s0, a1         # t2 = i * n
    slli    t2, t2, 2          # t2 = i * n * 4
    add     t3, a0, t2         # t3 = base address of row i
    slli    t4, s1, 2          # t4 = j * 4
    add     t5, t3, t4         # t5 = address of matrix[i][j]

    # Load vector of row elements starting at matrix[i][j]
    vle32.v v0, (t5)           # v0 <- matrix[i][j ... j+vl-1]

    # --- Compute index vector for loading column elements ---
    # For each k in 0..vl-1, we want the address for matrix[j+k][i]:
    # Address = a0 + [ (j+k)*n*4 + i*4 ]
    slli    t6, a1, 2          # t6 = n * 4 (stride for one row)
    slli    s7, s0, 2          # s7 = i * 4 (column offset)

    # Create a vector of indices: v2 = {0,1,...,vl-1}
    vid.v     v2                 # v2[i] = i for i=0...vl-1
    vadd.vx v2, v2, s1         # v2 = {j, j+1, ..., j+vl-1}
    vmul.vx v2, v2, t6         # v2 = {j*n*4, (j+1)*n*4, ...}
    vadd.vx v2, v2, s7         # v2 = {j*n*4 + i*4, (j+1)*n*4 + i*4, ...}

    # Use indexed load to fetch column elements: matrix[j...][i]
    vluxei32.v v1, (a0), v2       # v1 <- { matrix[j][i], matrix[j+1][i], ... }

    # --- Swap the two sets of elements -----------------------
    vse32.v v1, (t5)           # Store v1 into row i at matrix[i][j ...]
    vsuxei32.v v0, (a0), v2       # Store v0 into column i at matrix[j...][i]

    # Update inner loop index: j = j + vl (t0 holds vl)
    add     s8, s1, t0
    mv      s1, s8
    j       inner_loop

next_i:
    addi    s0, s0, 1          # i++
    j       outer_loop

done:
    ret
## END YOU CODE HERE

# Function: print
# Logs values from array in a0 into registers v1 for debugging and output.
# Inputs:
#   - a0: Base address of array
#   - a1: Size of array i.e. number of elements to log
# Clobbers: t0,t1, t2,t3 ft0, ft1.
printToLogVectorized:        
    addi sp, sp, -4
    sw a0, 0(sp)

    li t0, 0x123                 # Pattern for help in python script
    li t0, 0x456                 # Pattern for help in python script
    mv a1, a1                   # moving size to get it from log 
    mul a1, a1, a1              # sqaure matrix has n^2 elements 
	li t0, 0		                # load i = 0
    printloop:
        vsetvli t3, a1, e32           # Set VLEN based on a1
        slli t4, t3, 2                # Compute VLEN * 4 for address increment

        vle32.v v1, (a0)              # Load real[i] into v1
        add a0, a0, t4                # Increment pointer for real[] by VLEN * 4
        add t0, t0, t3                # Increment index

        bge t0, a1, endPrintLoop      # Exit loop if i >= size
        j printloop                   # Jump to start of loop
    endPrintLoop:
    li t0, 0x123                    # Pattern for help in python script
    li t0, 0x456                    # Pattern for help in python script
	
    lw a0, 0(sp)
    addi sp, sp, 4

	jr ra



# Function: _finish
# VeeR Related function which writes to to_host which stops the simulator
_finish:
    li x3, 0xd0580000
    addi x5, x0, 0xff
    sb x5, 0(x3)
    beq x0, x0, _finish

    .rept 100
        nop
    .endr


.data
## ALL DATA IS DEFINED HERE LIKE MATRIX, CONSTANTS ETC


## DATA DEFINE START
.equ MatrixSize, 5
matrix:
    .float -10.0, 13.0, 10.0, -3.0, 2.0
    .float 6.0, 15.0, 4.0, 13.0, 4.0
    .float 18.0, 2.0, 9.0, 8.0, -4.0
    .float 5.0, 4.0, 12.0, 17.0, 6.0
    .float -10.0, 7.0, 13.0, -3.0, 160.0
## DATA DEFINE END
size: .word MatrixSize