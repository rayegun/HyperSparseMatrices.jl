module HyperSparseMatrices
using SparseArrays

export HyperSparseMatrix, HyperSparseCSC, HyperSparseCSR
# This shouldn't be here. It should be in (imo) ArrayInterface.jl or a future AbstractSparse.jl
abstract type StorageOrder end
struct ColMajor <: StorageOrder end #colexicographic ordering
struct RowMajor <: StorageOrder end #lexicographic ordering

storageorder(::Array) = ColMajor()
storageorder(A::AbstractArray) = storageorder(parent(A))



struct HyperSparseMatrix{O, Tv, Tf, Ti<:Integer} <: AbstractSparseMatrix{Tv, Ti}
    # Comments here reflect column major ordering. 
    # To understand as DCSR simply swap references to columns with rows and vice versa.
    vlen::Int # m in CSC, n in CSR. This is the length of the vectors. In DCSC this is thus the number of rows.
    vdim::Int # n in CSC, m in CSR. This is the number of vectors being stored. in DCSC this is the number of columns.
    p::Vector{Ti} # The pointers into i/nzval. The row indices found in the k'th stored column are
    # found in idx[p[k], p[k+1]-1]
    h::Vector{Ti} # If column (row) j has stored entries, then j = h[k] for some k.
    # that is j is the k'th stored column.
    idx::Vector{Ti} # the stored row indices.
    v::Vector{Tv} # the coefficients of stored indices in the matrix
    fill::Tf
end

storageorder(::HyperSparseMatrix{O}) where {O} = O

const HyperSparseCSC{Tv, Tf, Ti} = HyperSparseMatrix{ColMajor(), Tv, Tf, Ti}
const HyperSparseCSR{Tv, Tf, Ti} = HyperSparseMatrix{RowMajor(), Tv, Tf, Ti}

Base.size(A::HyperSparseCSC) = (A.vlen, A.vdim)
Base.size(A::HyperSparseCSR) = (A.vdim, A.vlen)

SparseArrays.nnz(A::HyperSparseMatrix) = length(A.v)
SparseArrays.nonzeros(A::HyperSparseMatrix) = A.v

nvec(A::HyperSparseCSC) = size(A, 2)
nvec(A::HyperSparseCSR) = size(A, 1)

# indexing adapted from SS:GrB and SparseArrays:
function Base.getindex(A::HyperSparseMatrix{O, Tv, Tf, Ti}, row::Integer, col::Integer) where {O, Tv, Tf, Ti}
    # A HyperSparseMatrix stores sparse vectors.
    # Since we support both ColMajor and RowMajor i and j don't mean what they usually do.
    # Instead j refers to the vector
    # and i refers a scalar within that vector.
    if O === ColMajor() # For column major this is as it is for SparseMatrixCSC
        i = row
        j = col
    elseif O === RowMajor() # For row major, j is a row vector (which in this case (DCSC), may be compressed out entirely).
        i = col
        j = row
    end
    @boundscheck checkbounds(A, row, col) # the overloaded size should take care of doing this correctly.
    nnz(A) == o && return A.fill # no values, return fill.

    # A.h contains a list of vectors that were not compressed out.
    # We determine the idx in A.h of j (or the next greater idx if j ∉ A.h)
    idx = searchsortedfirst(A.h, j, 1, nvec(A), Forward) 

    # if this is true then we have a vector, and this just becomes the same as indexing in SparseMatrixCSC
    # If we don't, we just return the fill-in value.
    j == A.h[idx] || return A.fill 

    # this is all adapted directly from sparsematrix.jl ∈ SparseArrays.jl
    r1 = A.p[idx]
    r2 = A.p[idx + 1] - 1
    (r1 > r2) && return A.fill
    r1 = searchsortedfirst(A.idx, i, r1, r2, Forward)
    return ((r1 > r2) || (A.idx[r1] != i)) ? A.fill : A.v[r1]
end



#TODO:
# indexing
# printing (braille?)
# implement a good subset of the AbstractArray and AbstractSparseMatrix interface.
end
