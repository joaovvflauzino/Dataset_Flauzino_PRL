using DelimitedFiles, Combinatorics

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------
 orderN = 4                       # Size of the binary matrix (NxN)
 Nmatrices = 2^(orderN^2)        # Total number of possible binary matrices

# ------------------------------------------------------------
# FUNCTION: PermuteLines
# Permutes the rows of a matrix according to a given permutation.
# ------------------------------------------------------------
function permute_lines(M::AbstractMatrix{<:Integer}, perm::Vector{Int})
    return M[perm, :]                  # Return permuted matrix (no copy needed, Julia handles views efficiently)
end

# ------------------------------------------------------------
# FUNCTION: matrix_to_vector
# Converts a matrix to a column-major vector
# ------------------------------------------------------------
function matrix_to_vector(M::AbstractMatrix{<:Integer})
    return vec(M)                      # Flatten matrix in column-major order
end

# ------------------------------------------------------------
# FUNCTION: vector_to_matrix
# Converts a binary string or vector to a matrix
# ------------------------------------------------------------
function vector_to_matrix(binvec::AbstractString, N::Int)
    bits = parse.(Int, collect(binvec))       # Convert string to vector of bits
    M = reshape(bits, N, N)                   # Column-major reshape
    return M
end

# ------------------------------------------------------------
# FUNCTION: idx_label
# Computes a unique integer label for a binary NxN matrix
# using row-major order with (1,1) as the most significant bit.
# ------------------------------------------------------------
function idx_label(RM::AbstractMatrix{<:Integer})
    bits = vec(RM)'                            # Flatten matrix in row-major
    bitstring = join(bits)                     # Concatenate as string
    return parse(Int, bitstring; base=2)      # Convert binary string to integer
end

# ------------------------------------------------------------
# PREPARE LABEL STORAGE
# ------------------------------------------------------------
RMlabel = zeros(Int, Nmatrices)               # Initialize labels for all matrices
class = 1                                      # class counter

# ------------------------------------------------------------
# GENERATE ALL PERMUTATIONS OF ROWS
# Used for row/column swapping to identify equivalent matrices
# ------------------------------------------------------------
all_perms = collect(permutations(1:orderN))  # Symmetric class S_N

# ------------------------------------------------------------
# MAIN LOOP: iterate over all binary matrices
# Assign a unique class label to each equivalence class
# under row/column permutations and transpositions
# ------------------------------------------------------------
for i in 0:Nmatrices-1
    global class  
    if RMlabel[i+1] == 0                     # Skip already labeled matrices
        # Convert integer to binary string of length orderN^2
        binstr = bitstring(i)[end-orderN^2+1:end]
        RM = vector_to_matrix(binstr, orderN)

        # Iterate over all row permutations
        for perm_row in all_perms
            RM1 = permute_lines(RM, perm_row)

            # Iterate over all column permutations via transpose
            RMt = transpose(RM1)
            for perm_col in all_perms
                RM2 = permute_lines(RMt, perm_col)

                # Label the matrix and its transpose equivalence
                RMlabel[idx_label(RM2)+1] = class
                RMlabel[idx_label(transpose(RM2))+1] = class
            end
        end

        class += 1                             # Increment class counter
    end
end

# ------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------
writedlm("Labels_disorder_N4.txt", RMlabel)
