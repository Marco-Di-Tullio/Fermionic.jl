using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test basis(Op(4))[24] == 1.0
    @test basis(Op(6))[1] == 0.0
    @test cmcd(Op(4),1,2)[5,9] == -1.0
    @test cmcm(Op(4),4,2)[9,14] == 1.0
    @test cdcd(Op(4),2,3)[15,9] == 1.0
    @test vacuum(Op(4))[1] == 1.0
    @test basis_m(8,2)[2][4] == 10.0
    @test cdc(6,2,3,3)[6,6] == 1.0
    @test cdc(6,2,3,2)[14,15] == 1.0
    @test cdc(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[4,2] == 1
    @test ccd(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[2,4] == 1
    @test ccd(6,2,1,2)[8,12] == 1.0
    @test cdc(basis_m(4,2)[1],basis_m(4,2)[2],3,3)[5,5] == 1
    @test fixed_state([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],2)[6] == 1
    @test st(fixed_state(State([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],Op(4)),2))[6] == 1
    @test unfixed_state([0,0,0,0,1/sqrt(2),1/sqrt(2)],4,2)[13] == 1/sqrt(2)
    @test st(unfixed_state(State_fixed([0,0,0,0,0,1],Op(4),2)))[13] == 1
end
