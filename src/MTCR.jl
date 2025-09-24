module MTCR

using DocStringExtensions
using Libdl
using DifferentialEquations
using DiffEqBase: DiscreteCallback, set_proposed_dt!, u_modified!
using Printf

const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

include("data_conversion.jl")
include("fortran_wrapper.jl")
include("mtcr_config.jl")
include("mtcr_solver.jl")

export initialize_mtcr, finalize_mtcr
export MTCRConfig, MTCRResults
export solve_mtcr_0d, nitrogen_10ev_example

end
