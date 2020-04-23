
#SafeTestsets makes every test run in a separate enviromment
#so error can be found independently
using SafeTestsets
@safetestset "My f tests" begin include("my_f_tests.jl") end
