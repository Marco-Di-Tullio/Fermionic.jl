using Fermionic
using SparseArrays
using Test

@testset "Fermionic.jl" begin
    @test ssp(State(spzeros(16), Op(4))) == 0
    @test ssp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test ssp(State([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test ssp(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4))) == 1.0
    @test ssp(State_fixed(1/sqrt(2)*[1,0,0,0,0,1],Op_fixed(4,2))) == 1.0
    @test ssp(State_fixed(1/sqrt(2)*[1,0,0,0,0,im],Op_fixed(4,2))) == 1.0
    @test ssp(State_fixed(1/sqrt(2)*[1,0,0,0,0,1],Op_fixed(4,2))) == 1.0
    @test ssp(State_fixed(spzeros(6),Op_fixed(4,2))) == 0
    @test ssp(State_fixed(SparseVector(sparse([1,6],[1,1],[1/sqrt(2),im/sqrt(2)])),Op_fixed(4,2))) == 1.0
    @test sqsp(State(spzeros(16), Op(4))) == 0
    @test sqsp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1.0
    @test sqsp(State([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4))) == 0.8112781244591327
    @test sqsp(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4))) == 1
    @test majorization_sp(State([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4)),State([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4))) == 0
    @test majorization_qsp(State([0,0,0,im/2,0,0,0,0,0,0,0,0,sqrt(3)/2,0,0,0],Op(4)),State([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4))) == 1
    @test n_avg(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4))) == 2.0
    @test coef(State(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,5)[6,1] == round(1/sqrt(2),digits=5)
    @test coef(State_fixed(fixed_state(1/sqrt(2)*(ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),4),Op_fixed(6,4)),2,6)[13,10] == round(1/sqrt(2),digits=6)
    @test rhomd(State(1/sqrt(2)*(ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,5)[3,3] == 0.5
    @test rhomd(State(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,8)[3,3] == 0.5
    @test rhomd(State_fixed(fixed_state(1/sqrt(2)*(ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),4),Op_fixed(6,4)),2,7)[6,6] == 0.5
    @test rhomd(State_fixed(fixed_state(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),4),Op_fixed(6,4)),2,7)[1,1] == 0.5
    @test rhom(State_fixed(fixed_state(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),4),Op_fixed(6,4)),2,8)[1,15] == -0.5im
    @test rhom(State(1/sqrt(2)*(ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,7)[15,1] == 0.5
    @test rhom(State(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,7)[4,4] == 0.5
    @test rhom(State(1/sqrt(2)*(ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),Op(6)),2,4)[1,15] == 0.5
    @test rhom2(State_fixed(fixed_state(1/sqrt(2)*(im*ad(Op(6),1)*ad(Op(6),2)*ad(Op(6),3)*ad(Op(6),4)+ad(Op(6),3)*ad(Op(6),4)*ad(Op(6),5)*ad(Op(6),6))*vacuum(Op(6)),4),Op_fixed(6,4)))[1,15] == -0.5im
    @test trp(State(spzeros(16), Op(4)),[1])[1,1] == 0
    @test round(trp(State([0,0,0,1/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4)),[2,4])[2,2],digits=5) == 0.5
    @test round(trp(State([0,0,0,im/sqrt(2),0,0,0,0,0,0,0,0,1/sqrt(2),0,0,0],Op(4)),[1,4])[2,2],digits=5) == 0.5
    @test round(trp(State(1/sqrt(2)*(im*ad(Op(4),1)*ad(Op(4),2)+ad(Op(4),3)*ad(Op(4),4))*vacuum(Op(4)),Op(4)),[2,3])[2,2],digits=5) == 0.5
end
