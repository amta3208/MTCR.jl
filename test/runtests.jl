using MTCR: MTCR as mtcr
using Test
using Aqua

"""
Reset the MTCR state (Fortran + Julia) and initialize from scratch.

- `lib_path`: Filesystem path to the MTCR shared library
- `case_path`: Path to a case directory containing `input/prob_setup.inp`
- `config` (keyword, optional): If provided, a temporary case directory is created,
  input files are generated from this config there, and initialization is performed
  from that temporary directory (leaving `case_path` untouched).

Returns a NamedTuple with `(num_species, num_dimensions)`.
"""
function reset_and_init!(lib_path::AbstractString, case_path::AbstractString;
        config::Union{Nothing, mtcr.MTCRConfig} = nothing)
    # Best-effort cleanup (ignore errors if library is not loaded yet)
    try
        mtcr.finalize_api_wrapper()
    catch
        # ignore
    end
    try
        mtcr.close_mtcr_library()
    catch
        # ignore
    end

    # Fresh load + init
    mtcr.load_mtcr_library!(lib_path)
    # For test stability within one process, do not have MTCR finalize MPI
    try
        mtcr.set_api_finalize_mpi_wrapper(false)
    catch
        # ignore if symbol missing (older library)
    end

    if config === nothing
        return mtcr.initialize_api_wrapper(case_path = case_path)
    else
        # Use a temporary case directory so we don't touch the shared test_case
        mktempdir() do tmp
            # Generate inputs in the temp directory based on the provided config
            mtcr.generate_input_files(config, tmp)
            # Initialize from the temp directory
            dims = mtcr.initialize_api_wrapper(case_path = tmp)
            return dims
        end
    end
end

# Path to the MTCR shared library for tests
# TODO: Update shared library handling
const temp_mtcr_path = "/Users/amin/.julia/dev/MTCR/mtcr/source/libmtcr.so"

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

    # @testset "MTCR Solver" begin
    #     include("mtcr_solver.jl")
    # end
end
