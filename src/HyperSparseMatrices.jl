module HyperSparseMatrices
using SparseArrays
using StorageOrders
export HyperSparseMatrix, DCSCMatrix, DCSRMatrix


struct HyperSparseMatrix{O, Bi, Tv, Tf, Ti<:Integer} <: AbstractSparseMatrix{Tv, Ti}
    # Comments here reflect column major ordering. 
    # To understand as DCSR simply swap references to columns with rows and vice versa.
    vlen::Int # m in CSC, n in CSR. This is the length of the vectors. In DCSC this is thus the number of rows.
    vdim::Int # n in CSC, m in CSR. This is the number of vectors being stored. in DCSC this is the number of columns.

    p::Vector{Ti} # The pointers into i/nzval. The row indices found in the k'th stored column are
    # found in idx[p[k], p[k+1]-1]

    # If column (row) j has stored entries, then j = h[k] for some k.
    # j is the k'th stored column.
    h::Vector{Ti} 
    idx::Vector{Ti} # the stored row indices.
    v::Vector{Tv} # the coefficients of stored indices in the matrix
    fill::Tf # the fill value. A[i,j] == fill if (i,j) is compressed out.
end

StorageOrders.storageorder(::HyperSparseMatrix{O}) where {O} = O

const DCSCMatrix{Bi, Tv, Tf, Ti} = HyperSparseMatrix{ColMajor(), Bi, Tv, Tf, Ti}
const DCSRMatrix{Bi, Tv, Tf, Ti} = HyperSparseMatrix{RowMajor(), Bi, Tv, Tf, Ti}

DCSCMatrix{Bi}(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf} = DCSCMatrix{Tv, Tf, Ti}(vlen, vdim, p, h, idx, v, fill)
DCSRMatrix{Bi}(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf} = DCSRMatrix{Tv, Tf, Ti}(vlen, vdim, p, h, idx, v, fill)

DCSCMatrix(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf} = DCSCMatrix{1}(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf}
DCSRMatrix(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf} = DCSRMatrix{1}(vlen, vdim, p::Ti, h::Ti, idx::Ti, v::Tv, fill::Tf) where {Ti, Tv, Tf}

Base.size(A::DCSCMatrix) = (A.vlen, A.vdim)
Base.size(A::DCSRMatrix) = (A.vdim, A.vlen)

SparseArrays.nnz(A::HyperSparseMatrix) = length(A.v)
SparseArrays.nonzeros(A::HyperSparseMatrix) = A.v

nvec(A::DCSCMatrix) = size(A, 2)
nvec(A::DCSRMatrix) = size(A, 1)

# Bi setup taken from SparseMatrixCSR.jl
getoffset(S::HyperSparseMatrix{O, Bi}) where {O, Bi} = getoffset(Bi)
@inline getoffset(Bi::Integer) = 1-Bi

# indexing adapted from SS:GrB and SparseArrays:
function Base.getindex(A::HyperSparseMatrix{O, Bi, Tv, Tf, Ti}, row::Integer, col::Integer) where {O, Bi, Tv, Tf, Ti}
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
    o = getoffset(A)
    # A.h contains a list of vectors that were not compressed out.
    # We determine the idx in A.h of j (or the next greater idx if j ∉ A.h)
    io = i1 - o
    jo = j1 - o
    idx = searchsortedfirst(A.h, jo, 1, nvec(A), Forward) 

    # if this is true then we have a vector, and this just becomes the same as indexing in SparseMatrixCSC
    # If we don't, we just return the fill value.
    j == A.h[idx+o] || return A.fill 

    # this is all adapted directly from sparsematrix.jl ∈ SparseArrays.jl
    r1 = A.p[idx] + o
    r2 = A.p[idx + 1] - Bi
    (r1 > r2) && return A.fill
    r1 = searchsortedfirst(A.idx, io, r1, r2, Forward)
    return ((r1 > r2) || (A.idx[r1] != io)) ? A.fill : A.v[r1]
end



#TODO:
# COO => hypersparse
# hypersparse => COO
# tests
# indexing
# printing (braille?)
# implement a good subset of the AbstractArray and AbstractSparseMatrix interface.
end
