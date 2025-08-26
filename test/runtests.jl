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

    @testset "Fortran Wrapper" begin
        include("fortran_wrapper.jl")
    end

    @testset "MTCR Configuration" begin
        include("mtcr_config.jl")
    end
end
