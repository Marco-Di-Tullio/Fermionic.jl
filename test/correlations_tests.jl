using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ssp(State_sparse(spzeros(16), Op(4))) == 0
    @test ssp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test ssp(State_complex([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test ssp(State_sparse_complex(1/sqrt(2)*(im*cdm(Op(4),1)*cdm(Op(4),2)+cdm(Op(4),3)*cdm(Op(4),4))*vacuum(Op(4)),Op(4))) == 1.0
    @test sqsp(State_sparse(spzeros(16), Op(4))) == 0
    @test sqsp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test sqsp(State_complex([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4))) == 0.8112781244591327
    @test sqsp(State_sparse_complex(1/sqrt(2)*(im*cdm(Op(4),1)*cdm(Op(4),2)+cdm(Op(4),3)*cdm(Op(4),4))*vacuum(Op(4)),Op(4))) == 1
    @test majorization_sp(State_complex([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4)),State_complex([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4))) == 0
    @test majorization_qsp(State_complex([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4)),State_complex([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1
end
