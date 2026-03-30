# =============================================================================
# Code for quantifying disorder in data
#
# If you use this code in your research, please cite the following article
# where the methodology and theoretical formulation were introduced:
#
#   Flauzino, J. V. V., Prado, T. L., Marwan, N., Kurths, J., & Lopes, S. R.
#   “Quantifying Disorder in Data,” Phys. Rev. Lett. 135, 097401 (2025),
#   DOI: https://doi.org/10.1103/1y98-x33s
#
# This citation acknowledges the original development of the disorder
# quantifier based on recurrence analysis.
#
# This module implements an sliding-window algorithm to count all
# 4×4 binary matrices derived from recurrence microstates of a (multivariate) time series
# and use it as input to quantify disorder in time series.
#
# Go to the end of the code to see an example of its use.
# =============================================================================

module FastMicrostates4

using Base.Threads

# ============================================================
# This module computes the full set of 4x4 binary microstates
# extracted from recurrence matrices obtained from a time series.
#
# Convention:
# - 4x4 binary matrices
# - Row-major ordering
# - Element (1,1) is the most significant bit (MSB)
# ============================================================


# ------------------------------------------------------------
# Computes the integer label of a 4x4 binary matrix using
# EXACTLY the same convention as the external label generator:
#
# - Row-major order
# - (1,1) is the most significant bit
# - (4,4) is the least significant bit
# ------------------------------------------------------------
@inline function idxlabel_fast(RM::AbstractMatrix)
    size(RM, 1) == 4 && size(RM, 2) == 4 ||
        throw(ArgumentError("RM must be 4x4"))

    code = UInt16(0)

    @inbounds for i in 1:4
        for j in 1:4
            # Shift left and insert current bit
            code = (code << 1) | (iszero(RM[i, j]) ? UInt16(0) : UInt16(1))
        end
    end

    return Int(code)
end


# ------------------------------------------------------------
# Computes the full pairwise squared Euclidean distance matrix.
#
# Input:
#   X : K × d matrix (K time points, d dimensions)
#
# Output:
#   D2 : K × K symmetric matrix of squared distances
#
# The computation is parallelized over rows.
# ------------------------------------------------------------
function pairwise_sqeuclidean(X::AbstractMatrix{<:Real})
    Xf = Matrix{Float64}(X)
    K, d = size(Xf)
    D2 = Matrix{Float64}(undef, K, K)

    @threads for i in 1:K
        D2[i, i] = 0.0
        @inbounds for j in i+1:K
            s = 0.0
            @simd for k in 1:d
                δ = Xf[i, k] - Xf[j, k]
                s += δ * δ
            end
            D2[i, j] = s
            D2[j, i] = s
        end
    end

    return D2
end


# ------------------------------------------------------------
# Builds the 4-bit horizontal codes for one row of the
# thresholded recurrence matrix.
#
# Bit ordering inside the 4-bit block:
#   column 1 -> bit 3
#   column 2 -> bit 2
#   column 3 -> bit 1
#   column 4 -> bit 0
#
# This matches the row-major MSB convention.
# ------------------------------------------------------------
@inline function fill_rowcodes!(
    rowcodes::Vector{UInt16},
    D2::AbstractMatrix{T},
    i::Int,
    τ::T
) where {T<:AbstractFloat}

    W = length(rowcodes)

    @inbounds begin
        code = UInt16(0)

        # Build first 4-bit window explicitly
        if D2[i, 1] <= τ; code |= 0x0008; end
        if D2[i, 2] <= τ; code |= 0x0004; end
        if D2[i, 3] <= τ; code |= 0x0002; end
        if D2[i, 4] <= τ; code |= 0x0001; end

        rowcodes[1] = code

        # Sliding window horizontally
        for j in 2:W
            newbit = D2[i, j + 3] <= τ ? UInt16(1) : UInt16(0)

            # Keep lower 3 bits, shift left, insert new bit
            code = ((code & 0x0007) << 1) | newbit

            rowcodes[j] = code
        end
    end

    return nothing
end


# ------------------------------------------------------------
# Accumulates all 4x4 microstates for one vertical window.
#
# The final 16-bit code is constructed row-major:
#
#   r1 -> most significant 4 bits
#   r4 -> least significant 4 bits
#
# Omega[label+1] is incremented.
# ------------------------------------------------------------
@inline function accumulate_windows!(
    Omega::AbstractVector{UInt32},
    r1::Vector{UInt16},
    r2::Vector{UInt16},
    r3::Vector{UInt16},
    r4::Vector{UInt16}
)
    W = length(r1)

    @inbounds for j in 1:W
        code = (UInt16(r1[j]) << 12) |
               (UInt16(r2[j]) << 8)  |
               (UInt16(r3[j]) << 4)  |
                UInt16(r4[j])

        Omega[Int(code) + 1] += 1
    end

    return nothing
end


# ------------------------------------------------------------
# Counts microstates for a single epsilon value.
#
# Uses a rolling buffer of 4 rows to avoid recomputation.
# ------------------------------------------------------------
function count_one_eps_from_sqdist!(
    Omega_col::AbstractVector{UInt32},
    D2::AbstractMatrix{T},
    eps::T
) where {T<:AbstractFloat}

    K = size(D2, 1)
    K >= 4 || throw(ArgumentError("Distance matrix must have at least 4 points."))

    W = K - 3
    τ = eps * eps   # squared threshold

    # Four rolling buffers for horizontal codes
    bufs = [Vector{UInt16}(undef, W) for _ in 1:4]

    # Initialize first four rows
    fill_rowcodes!(bufs[1], D2, 1, τ)
    fill_rowcodes!(bufs[2], D2, 2, τ)
    fill_rowcodes!(bufs[3], D2, 3, τ)
    fill_rowcodes!(bufs[4], D2, 4, τ)

    s1, s2, s3, s4 = 1, 2, 3, 4

    @inbounds for top in 1:W
        accumulate_windows!(Omega_col, bufs[s1], bufs[s2], bufs[s3], bufs[s4])

        if top < W
            # Reuse buffer of row leaving window
            fill_rowcodes!(bufs[s1], D2, top + 4, τ)

            # Rotate buffers
            s1, s2, s3, s4 = s2, s3, s4, s1
        end
    end

    return nothing
end


# ------------------------------------------------------------
# Main user interface.
#
# Receives time series X (K × d).
# Automatically generates epsilon values in the range:
#   [1e-6, 0.5 * d_max]
#
# Returns:
#   Omega[label+1, epsilon_index]
# ------------------------------------------------------------
function count_microstates4(
    X::AbstractMatrix{<:Real};
    neps::Int = 40
)
    D2 = pairwise_sqeuclidean(X)

    # Maximum actual distance
    dmax = sqrt(maximum(D2))

    # Epsilon vector
    epsilons = collect(range(1e-6, 0.5 * dmax, length = neps))

    Omega = count_microstates4_from_sqdist(D2, epsilons)

    return Omega
end


# ------------------------------------------------------------
# Same as above, but receives squared distance matrix directly.
# ------------------------------------------------------------
function count_microstates4_from_sqdist(
    D2::AbstractMatrix{<:Real},
    epsilons::AbstractVector{<:Real}
)
    D2f = Matrix{Float64}(D2)
    E = length(epsilons)

    # 2^(4*4) = 65536 total microstates
    Omega = zeros(UInt32, 65536, E)

    # Parallel over epsilon values
    @threads for e in 1:E
        eps = Float64(epsilons[e])
        count_one_eps_from_sqdist!(view(Omega, :, e), D2f, eps)
    end

    return Omega
end

end # module


# ======================================================================
# Disorder quantification using precomputed recurrence microstate labels
# ======================================================================
function QuantifyDisorder(X::AbstractMatrix{<:Real}) where T<:Real

    dimension = length(X[1, :])

    # Compute microstate counts
    Omega = FastMicrostates4.count_microstates4(X)

    # Load equivalence classes of 4x4 matrices computed using Eq.8-Flauzino-PRL-2025
    RMlabel = parse.(Int, readlines("Labels_disorder_N4.txt"))

    labels = @view RMlabel[:, 1]
    nclasses = labels[end]
    leneps = size(Omega, 2)
    Nmatrices = size(Omega, 1)
    
    # Build index list of each equivalence class
    class_idxs = [Int[] for _ in 1:nclasses]
    for i in 1:Nmatrices
        g = labels[i]
        push!(class_idxs[g], i)
    end

    Entropy = zeros(Float64, leneps)

    # Loop over equivalence classes
    for g in 1:nclasses
        idxs = class_idxs[g]
        sizeclass = length(idxs)
        sizeclass > 1 || continue

        # Normalization factor per class
        norm = vec(sum(Omega[idxs, :], dims=1))
        logsize = log(sizeclass)

        for j in 1:leneps
            normj = norm[j]
            normj == 0.0 && continue

            s = 0.0
            for i in idxs
                val = Omega[i, j]
                val == 0.0 && continue
                p = val / normj
                s -= p * log(p)
            end

            Entropy[j] += s / logsize
        end
    end

    #Normalization constant
    A = 190
    if dimension == 1
        A = 145
    end

    return Entropy ./ A
end


# ============================================================
# EXAMPLE USAGE
# ============================================================
X = randn(3000, 3)
Disorder = maximum(QuantifyDisorder(X))
