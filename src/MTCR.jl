module MTCR

using DocStringExtensions
using DelimitedFiles: readdlm, writedlm

const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

include("data_conversion.jl")

end
