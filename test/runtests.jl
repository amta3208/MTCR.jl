using MTCR
using Test
using Aqua

@testset "MTCR.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MTCR; ambiguities = false,)
    end
    # Write your tests here.
end
