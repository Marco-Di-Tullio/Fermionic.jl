using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ssp(State([0 for i in 1:16], Op(4))) == 0
    #in order to execute these tests, enter the pkg rpl with ]
    #and execute "test PackageName"
end
