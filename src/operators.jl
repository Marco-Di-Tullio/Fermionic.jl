#operators constructors

struct Op
    dim::Int
    cmtot::SparseMatrixCSC{Float64,Int64}
    cdtot::SparseMatrixCSC{Float64,Int64}
    le::Int
    basis::SparseMatrixCSC{Float64,Int64}
    Op(dim) = new(dim, operators(dim)...)
end

#after being defined, one can
#retrieve the dimension by typping
# dim(op4)
dim(o::Op) = o.dim
cmtot(o::Op) = o.cmtot
cdtot(o::Op) = o.cdtot
le(o::Op) = o.le
basis(o::Op) = o.basis

#These are the fermionic operators we wanted
#By calling cm(op, 1) we get the sparse matrix
#corresponding to the destruction of the
#first fermionic mode
cm(o::Op, i::Int) = cmtot(o)[1:le(o),((i-1)*le(o)+1):i*le(o)]
cdm(o::Op, i::Int) = cdtot(o)[((i-1)*le(o)+1):i*le(o),1:le(o)]
cdcm(o::Op, i::Int, j::Int) = cdm(o,i)*cm(o,j)
cmcd(o::Op, i::Int, j::Int) = cm(o,i)*cdm(o,j)
cmcm(o::Op, i::Int, j::Int) = cm(o,i)*cm(o,j)
cdcd(o::Op, i::Int, j::Int) = cdm(o,i)*cdm(o,j)

function vacuum(o::Op)
    l = 2^dim(o)
    vac = spzeros(l)
    vac[1] = 1.0
    return vac
end
