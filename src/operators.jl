#operators constructors

struct Op
    dim::Int
    cmtot::SparseMatrixCSC{Float64,Int64}
    cdtot::SparseMatrixCSC{Float64,Int64}
    le::Int
    basis::SparseMatrixCSC{Float64,Int64}
    Op(dim) = new(dim, operators(dim)...)
end

# after being defined, one can
# retrieve the dimension by typping
# dim(o)
dim(o::Op) = o.dim
cmtot(o::Op) = o.cmtot
cdtot(o::Op) = o.cdtot
le(o::Op) = o.le
basis(o::Op) = o.basis

#These are the fermionic operators we wanted
#By calling cm(op, 1) we get the sparse matrix
#corresponding to the destruction of the
#first fermionic mode
a(o::Op, i::Int) = cmtot(o)[1:le(o),((i-1)*le(o)+1):i*le(o)]
ad(o::Op, i::Int) = cdtot(o)[((i-1)*le(o)+1):i*le(o),1:le(o)]
ada(o::Op, i::Int, j::Int) = ad(o,i)*a(o,j)
aad(o::Op, i::Int, j::Int) = a(o,i)*ad(o,j)
aa(o::Op, i::Int, j::Int) = a(o,i)*a(o,j)
adad(o::Op, i::Int, j::Int) = ad(o,i)*ad(o,j)

function vacuum(o::Op)
    l = 2^dim(o)
    vac = spzeros(l)
    vac[1] = 1.0
    return vac
end
