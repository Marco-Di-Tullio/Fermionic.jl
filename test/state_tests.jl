using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test rhoqsp(State(spzeros(16), Op(4)))[1,1] == 0
    @test rhoqsp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4)))[2,2] == 0.5
    @test rhoqsp(State([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4)))[1,2] == 0
    @test rhoqsp(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4)))[6,6] == 0.5
    @test rhoqsp(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4)),4)[6,6] == 0.5
    @test non_zero([0,0,0,1,0.5], 2)[2][1] == 1.0
    @test length(non_zero([0,0,0,1,0.00005], 4)[1]) == 1
    #in order to execute these tests, enter the pkg rpl with ]
    #and execute "test PackageName"
end
