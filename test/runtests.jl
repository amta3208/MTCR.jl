using MTCR: MTCR as mtcr
using Test
using Aqua

@testset "MTCR.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(mtcr; ambiguities = false)
    end

    @testset "Data Conversion" begin
        include("data_conversion.jl")
    end
end
