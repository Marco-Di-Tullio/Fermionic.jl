struct State2{T}
    st::Array{T,1}
    ope::Op
end

struct State_sparse2{T}
    st::SparseVector{T,Int64}
    ope::Op
end


st(s::State2) = s.st
st(s::State_sparse2) = s.st

typ(s::State2) = eltype(s.st)
