using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test dim(Op(8)) == 8
    @test basis(Op(4))[24] == 1.0
    @test basis(Op(6))[1] == 0.0
    @test aad(Op(4),1,2)[5,9] == -1.0
    @test ada(Op(4),3,2)[3,5] == 1
    @test aa(Op(4),4,2)[9,14] == 1.0
    @test adad(Op(4),2,3)[15,9] == 1.0
    @test vacuum(Op(4))[1] == 1.0
    @test basis_m(8,2)[2][4] == 10.0
    @test ada(6,2,3,3)[6,6] == 1.0
    @test ada(6,2,3,2)[14,15] == 1.0
    @test ada(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[4,2] == 1
    @test aad(basis_m(4,2)[1],basis_m(4,2)[2],1,2)[2,4] == -1
    @test aad(6,2,1,2)[8,12] == -1.0
    @test ada(basis_m(4,2)[1],basis_m(4,2)[2],3,3)[5,5] == 1
    @test fixed(ad(Op(6),1)*a(Op(6),1),2)[11,11] == 1.0
    @test fixed_state([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],2)[6] == 1
    @test unfixed_state([0,0,0,0,1/sqrt(2),1/sqrt(2)],4,2)[13] == 1/sqrt(2)
    @test st(unfixed_state(State_fixed([0,0,0,0,0,1],Op_fixed(4,2))))[13] == 1
    @test ada(Op_fixed(4,2),1,2)[4,2] == 1
    @test aad(Op_fixed(4,2),1,2)[2,4] == 1
    @test dim(Op_fixed(8,2)) == 8
    @test nume(Op_fixed(8,2)) == 2
    @test basis(Op_fixed(4,2))[4,1] == 1
    @test ad(6,3,2)[13,14] == -1
    @test ad(6,3,7)[29,19] == -1
    @test a(6,3,2)[12,18] == -1
    @test a(4,3,5)[1] == 0
    @test ad(6,2,[1,2])[10,1] == -1
    @test a(6,3,[1,2])[2,18] == 1
    @test af(4,2)[1] == a(4,2,1)
    @test adf(4,2)[2] == ad(4,2,2)
    #@test Op_semifixed(6).d == 6
    #@test ad(Op_semifixed(6),4,2)[3,6] == -1
    #@test a(Op_semifixed(6),4,2)[13,12] == -1
    #@test a(Op_semifixed(6),4,7)[1] == 0
    #@test ad(Op_semifixed(6),2,[1,2])[10,1] == -1
    #@test a(Op_semifixed(6),4,[1,2])[3,12] == 1
end
