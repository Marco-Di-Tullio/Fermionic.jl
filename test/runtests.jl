
#SafeTestsets makes every test run in a separate enviromment
#so error can be found independently
using SafeTestsets

@safetestset "Operator tests" begin include("operator_tests.jl") end
@safetestset "State tests" begin include("state_tests.jl") end
