module MTCR

using DocStringExtensions
using Libdl

const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

include("data_conversion.jl")
include("fortran_wrapper.jl")
include("mtcr_config.jl")

end
