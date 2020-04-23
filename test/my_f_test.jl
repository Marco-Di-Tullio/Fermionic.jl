using Fermionic
using Test


@testset "Fermionic.jl" begin
    @test my_f(3,1) == 7
    @test my_f(2,1) == 5
    #in order to execute these tests, enter the pkg rpl with ]
    #and execute "test PackageName"
end
