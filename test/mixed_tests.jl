using Fermionic
using SparseArrays
using Test

@testset "mixed.jl" begin
    @test eigensp_mixed([0.25,0.75],[State([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4)),State([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],Op(4))])[1] == 0.875
end
