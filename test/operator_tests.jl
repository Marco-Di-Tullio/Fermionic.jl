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
end
