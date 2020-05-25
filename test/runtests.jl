
#SafeTestsets makes every test run in a separate enviromment
#so error can be found independently
using SafeTestsets

@safetestset "Operator tests" begin include("operator_tests.jl") end
@safetestset "State tests" begin include("state_tests.jl") end
@safetestset "Logic gates tests" begin include("logic_gates_tests.jl") end
@safetestset "Correlations tests" begin include("correlations_tests.jl") end
@safetestset "Mixed tests" begin include("mixed_tests.jl") end
