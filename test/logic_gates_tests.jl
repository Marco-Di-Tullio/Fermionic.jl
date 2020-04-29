using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ucnot(Op(2),1,2)[4,3] == 1.0
    @test ucnot(Op(6),4,2)[48,64] == 1.0
    @test ucnot(Op(3),2,3)[1,2] == 0.0
end
