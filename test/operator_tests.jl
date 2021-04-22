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
    @test cdc(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[4,2]==1
    @test ccd(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[2,4]==1
end
