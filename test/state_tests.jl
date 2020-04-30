using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ssp(State_sparse(spzeros(16), Op(4))) == 0
    @test ssp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    #in order to execute these tests, enter the pkg rpl with ]
    #and execute "test PackageName"
end
