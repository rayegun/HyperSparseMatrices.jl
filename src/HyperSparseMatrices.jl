module HyperSparseMatrices
using SparseArrays

export HyperSparseMatrix, HyperSparseCSC, HyperSparseCSR
# This shouldn't be here. It should be in (imo) ArrayInterface.jl or a future AbstractSparse.jl
abstract type StorageOrder end
struct ColMajor <: StorageOrder end #colexicographic ordering
struct RowMajor <: StorageOrder end #lexicographic ordering

storageorder(::Array) = ColMajor()
storageorder(A::AbstractArray) = storageorder(parent(A))



struct HyperSparseMatrix{O, Tv, Ti<:Integer} <: AbstractSparseMatrix{Tv, Ti}
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
end

const HyperSparseCSC{Tv, Ti} = HyperSparseMatrix{ColMajor(), Tv, Ti}
const HyperSparseCSR{Tv, Ti} = HyperSparseMatrix{RowMajor(), Tv, Ti}

Base.size(A::HyperSparseCSC) = (A.vlen, A.vdim)
Base.size(A::HyperSparseCSR) = (A.vdim, A.vlen)

#TODO:
# indexing
# printing (braille?)
# implement a good subset of the AbstractArray and AbstractSparseMatrix interface.
end
