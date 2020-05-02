using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ucnot(Op(2),1,2)[4,3] == 1.0
    @test ucnot(Op(6),4,2)[48,64] == 1.0
    @test ucnot(Op(3),2,3)[1,2] == 0.0
    @test hadamard(Op(3),2,1)[2,2] == 1.0
    @test hadamard(Op(3),3,1)[5,2] == 1/sqrt(2)
    @test phase(Op(3),2,pi/4)[3,3] == round(1/sqrt(2)+im/sqrt(2),digits=15)
    @test phase(Op(6),2,pi/4)[1,3] == 0.0
    @test swap(Op(3),1,2)[3,5] == 1.0
    @test swap(Op(5),1,2)[5,6] == 0.0
    @test sigma_x(Op(3),2)[1,3] == 1.0
    @test sigma_x(Op(5),2)[1,1] == 0.0
    @test sigma_y(Op(2),2)[2,1] == 0.0 + 1.0*im
    @test sigma_y(Op(4),3)[2,4] == -1.0*im
    @test sigma_z(Op(4),3)[2,4] == 0.0
    @test sigma_z(Op(3),3)[3,3] == 1.0
end
