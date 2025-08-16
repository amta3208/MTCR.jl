"""
# Fortran Wrapper Module

This module provides low-level Julia bindings to the MTCR Fortran library using `ccall`.
It handles direct interfacing with the Fortran API functions defined in `interface.f90`.
"""

# MTCR library state
const MTCR_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const MTCR_LIB_PATH = Ref{String}("")
const MTCR_INITIALIZED = Ref{Bool}(false)

"""
$(SIGNATURES)

Set the path to the MTCR shared library and load it.

# Arguments
- `path::String`: Path to the MTCR shared library file

# Throws
- `ErrorException`: If the library cannot be loaded
"""
function load_mtcr_library!(path::String)
    # Close existing handle if open
    if MTCR_HANDLE[] != C_NULL
        Libdl.dlclose(MTCR_HANDLE[])
        MTCR_HANDLE[] = C_NULL
    end

    # Validate path exists
    if !isfile(path)
        error("MTCR library file not found: $path")
    end

    # Open new library
    try
        MTCR_HANDLE[] = Libdl.dlopen(path)
        MTCR_LIB_PATH[] = path  # Store the path for ccall usage
    catch e
        error("Failed to load MTCR library from $path: $(e.msg)")
    end
end

"""
$(SIGNATURES)

Get the handle to the loaded MTCR library.

# Returns
- `Ptr{Cvoid}`: Handle to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_mtcr_handle()
    if MTCR_HANDLE[] == C_NULL
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end
    return MTCR_HANDLE[]
end

"""
$(SIGNATURES)

Get the path to the loaded MTCR library.

# Returns
- `String`: Path to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_mtcr_lib_path()
    if MTCR_HANDLE[] == C_NULL
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end
    return MTCR_LIB_PATH[]
end

"""
$(SIGNATURES)

Check if the MTCR library is currently loaded.

# Returns
- `Bool`: True if library is loaded, false otherwise
"""
function is_mtcr_loaded()
    return MTCR_HANDLE[] != C_NULL
end

"""
$(SIGNATURES)

Close the MTCR library and free resources.
"""
function close_mtcr_library()
    if MTCR_HANDLE[] != C_NULL
        Libdl.dlclose(MTCR_HANDLE[])
        MTCR_HANDLE[] = C_NULL
        MTCR_LIB_PATH[] = ""  # Clear the path
    end
end

"""
$(SIGNATURES)

Initialize the MTCR API system.

# Arguments
- `case_path::String`: Path to directory containing input/ subdirectory (default: current directory)

# Returns
- `NamedTuple`: Contains `num_species` and `num_dimensions` as determined by MTCR from input files

# Throws
- `ErrorException`: If case_path doesn't exist, input file is missing, or Fortran call fails
"""
function initialize_api_wrapper(; case_path::String = pwd())
    if MTCR_INITIALIZED[]
        @warn "MTCR already initialized in this Julia session - skipping"
        return nothing
    end

    # Validate inputs
    if !isdir(case_path)
        error("Case path does not exist: $case_path")
    end

    input_file = joinpath(case_path, "input", "prob_setup.inp")
    if !isfile(input_file)
        error("Required input file not found: $input_file")
    end

    # Ensure output directory structure exists
    output_dir = joinpath(case_path, "output")
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Ensure required output subdirectories exist
    sources_dir = joinpath(output_dir, "sources")
    states_dir = joinpath(output_dir, "states")

    if !isdir(sources_dir)
        mkpath(sources_dir)
    end

    if !isdir(states_dir)
        mkpath(states_dir)
    end

    # Store current directory and change to case path
    original_dir = pwd()

    try
        cd(case_path)

        # Create references for the Fortran call (as output parameters)
        # Fortran interface expects (num_species, num_dimensions) order
        num_species_ref = Ref{Int32}(0)
        num_dimensions_ref = Ref{Int32}(0)

        # Call Fortran function with correct parameter order
        ccall((:initialize_api, get_mtcr_lib_path()), Cvoid,
            (Ref{Int32}, Ref{Int32}),
            num_species_ref, num_dimensions_ref)

        # If we get here, the call succeeded
        MTCR_INITIALIZED[] = true

        # Return the values determined by Fortran
        return (num_species = num_species_ref[], num_dimensions = num_dimensions_ref[])

    catch e
        # Re-throw with more context about the failure
        if isa(e, Base.SystemError) || isa(e, ErrorException)
            error("MTCR initialization failed in directory $case_path: $(e.msg)")
        else
            rethrow(e)
        end
    finally
        # Always restore original directory
        cd(original_dir)
    end
end

"""
$(SIGNATURES)

Finalize the MTCR API system and clean up resources.
"""
function finalize_api_wrapper()
    if !MTCR_INITIALIZED[]
        @warn "MTCR not initialized - nothing to finalize"
        return nothing
    end

    ccall((:finalize_api, get_mtcr_lib_path()), Cvoid, ())

    MTCR_INITIALIZED[] = false
    return nothing
end

"""
$(SIGNATURES)

Get the maximum number of species supported by MTCR.

# Returns
- `Int32`: Maximum number of species

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_max_number_of_species_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    return ccall((:get_max_number_of_species, get_mtcr_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the maximum number of electronic states per atomic species supported by MTCR.

# Returns
- `Int32`: Maximum number of atomic electronic states

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_max_number_of_atomic_electronic_states_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    return ccall(
        (:get_max_number_of_atomic_electronic_states, get_mtcr_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the maximum number of electronic states per molecular species supported by MTCR.

# Returns
- `Int32`: Maximum number of molecular electronic states

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_max_number_of_molecular_electronic_states_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    return ccall(
        (:get_max_number_of_molecular_electronic_states, get_mtcr_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get species names from MTCR.

# Returns
- `Vector{String}`: Array of species names

# Throws
- `ErrorException`: If MTCR library is not loaded
"""
function get_species_names_wrapper()
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    # Get max species count first
    max_species = get_max_number_of_species_wrapper()

    # Fortran character length (from parameters module: nmlen)
    name_length = 32  # This should match nmlen in Fortran

    # Allocate buffer for species names (Fortran returns character array)
    names_buffer = zeros(UInt8, name_length * max_species)

    # Call Fortran subroutine
    ccall((:get_species_names, get_mtcr_lib_path()), Cvoid,
        (Ptr{UInt8},), names_buffer)

    # Convert buffer to Julia strings
    species_names = String[]

    # The Fortran subroutine packs names sequentially with null terminators
    # Look for null-terminated strings
    i = 1
    while i <= length(names_buffer)
        # Find the next null terminator
        null_idx = findfirst(==(0), names_buffer[i:end])
        if null_idx === nothing
            break  # No more null terminators
        end

        # Extract the name bytes (excluding null terminator)
        name_end = i + null_idx - 2
        if name_end >= i
            name_bytes = names_buffer[i:name_end]
            name = String(name_bytes) |> strip
            if !isempty(name)
                push!(species_names, name)
            end
        end

        i += null_idx
        # Skip any additional null padding to get to next name
        while i <= length(names_buffer) && names_buffer[i] == 0
            i += 1
        end
    end

    return species_names
end

"""
$(SIGNATURES)

Calculate nonequilibrium source terms.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_etot::Float64`: Total energy density
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- Tuple of derivative arrays corresponding to input arrays

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
"""
function calculate_sources_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    # Prepare output arrays
    drho_sp = similar(rho_sp)
    drho_etot = Ref{Float64}(0.0)

    # Handle optional output arguments
    drho_ex = rho_ex !== nothing ? similar(rho_ex) : nothing
    drho_vx = rho_vx !== nothing ? similar(rho_vx) : nothing
    drho_erot = rho_erot !== nothing ? Ref{Float64}(0.0) : nothing
    drho_eeex = rho_eeex !== nothing ? Ref{Float64}(0.0) : nothing
    drho_evib = rho_evib !== nothing ? Ref{Float64}(0.0) : nothing

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_nonequilibrium_sources, get_mtcr_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # rho_u (optional)
                Ptr{Float64},                                    # rho_v (optional)
                Ptr{Float64},                                    # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64},                                    # rho_evib (optional)
                Ptr{Float64},                                    # drho_sp
                Ptr{Float64},                                    # drho_ex (optional)
                Ptr{Float64},                                    # drho_vx (optional)
                Ref{Float64},                                    # drho_etot
                Ptr{Float64},                                    # drho_erot (optional)
                Ptr{Float64},                                    # drho_eeex (optional)
                Ptr{Float64}),                                   # drho_evib (optional)
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            drho_sp,
            drho_ex !== nothing ? drho_ex : C_NULL,
            drho_vx !== nothing ? drho_vx : C_NULL,
            drho_etot,
            drho_erot !== nothing ? drho_erot : C_NULL,
            drho_eeex !== nothing ? drho_eeex : C_NULL,
            drho_evib !== nothing ? drho_evib : C_NULL)
    catch e
        error("Failed to calculate source terms: $(e)")
    end

    return (drho_sp = drho_sp,
        drho_etot = drho_etot[],
        drho_ex = drho_ex,
        drho_vx = drho_vx,
        drho_erot = drho_erot !== nothing ? drho_erot[] : nothing,
        drho_eeex = drho_eeex !== nothing ? drho_eeex[] : nothing,
        drho_evib = drho_evib !== nothing ? drho_evib[] : nothing)
end

"""
$(SIGNATURES)

Calculate temperatures from thermodynamic state.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities
- `rho_etot::Float64`: Total energy density
- Additional optional energy components

# Returns
- Named tuple with temperatures (tt, trot, teex, tvib, tex, tvx)

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
"""
function calculate_temperatures_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    # Prepare output variables
    tt = Ref{Float64}(0.0)
    trot = Ref{Float64}(0.0)
    teex = Ref{Float64}(0.0)
    tvib = Ref{Float64}(0.0)

    # Prepare arrays for species-specific temperatures
    max_species = get_max_number_of_species_wrapper()
    tex = zeros(Float64, max_species)
    # Get max molecular electronic states for tvx array
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    tvx = zeros(Float64, max_molecular_electronic_states, max_species)

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_temperatures, get_mtcr_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # rho_u (optional)
                Ptr{Float64},                                    # rho_v (optional)
                Ptr{Float64},                                    # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64},                                    # rho_evib (optional)
                Ref{Float64},                                    # tt
                Ref{Float64},                                    # trot
                Ref{Float64},                                    # teex
                Ref{Float64},                                    # tvib
                Ptr{Float64},                                    # tex
                Ptr{Float64}),                                   # tvx
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            tt,
            trot,
            teex,
            tvib,
            tex,
            tvx)
    catch e
        error("Failed to calculate temperatures: $(e)")
    end

    return (tt = tt[], trot = trot[], teex = teex[], tvib = tvib[],
        tex = tex, tvx = tvx)
end

"""
$(SIGNATURES)

Calculate total energy from state variables.

# Arguments
- `tt::Float64`: Translational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `u::Float64`: x-velocity component (optional)
- `v::Float64`: y-velocity component (optional)
- `w::Float64`: z-velocity component (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- `Float64`: Total energy density

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
"""
function calculate_total_energy_wrapper(tt::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        u::Union{Float64, Nothing} = nothing,
        v::Union{Float64, Nothing} = nothing,
        w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    rho_etot = Ref{Float64}(0.0)

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_total_energy, get_mtcr_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_etot (output)
                Ref{Float64},                                    # tt
                Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # rho_vx (optional)
                Ptr{Float64},                                    # u (optional)
                Ptr{Float64},                                    # v (optional)
                Ptr{Float64},                                    # w (optional)
                Ptr{Float64},                                    # rho_erot (optional)
                Ptr{Float64},                                    # rho_eeex (optional)
                Ptr{Float64}),                                   # rho_evib (optional)
            rho_etot,
            tt,
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            rho_vx !== nothing ? rho_vx : C_NULL,
            u !== nothing ? Ref{Float64}(u) : C_NULL,
            v !== nothing ? Ref{Float64}(v) : C_NULL,
            w !== nothing ? Ref{Float64}(w) : C_NULL,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL)
    catch e
        error("Failed to calculate total energy: $(e)")
    end

    return rho_etot[]
end

"""
$(SIGNATURES)

Calculate vibrational energy from state variables.

# Arguments
- `tvib::Float64`: Vibrational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `tex::Vector{Float64}`: Electronic temperatures per species (optional)
- `teex::Float64`: Electron-electronic temperature (optional)

# Returns
- `Float64`: Vibrational energy density

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
"""
function calculate_vibrational_energy_wrapper(tvib::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        tex::Union{Vector{Float64}, Nothing} = nothing,
        teex::Union{Float64, Nothing} = nothing)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    if tvib <= 0.0
        error("Invalid vibrational temperature: $tvib")
    end

    rho_evib = Ref{Float64}(0.0)

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_vibrational_energy, get_mtcr_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_evib (output)
                Ref{Float64},                                    # tvib
                Ptr{Float64},                                    # rho_sp
                Ptr{Float64},                                    # rho_ex (optional)
                Ptr{Float64},                                    # tex (optional)
                Ptr{Float64}),                                   # teex (optional)
            rho_evib,
            tvib,
            rho_sp,
            rho_ex !== nothing ? rho_ex : C_NULL,
            tex !== nothing ? tex : C_NULL,
            teex !== nothing ? Ref{Float64}(teex) : C_NULL)
    catch e
        error("Failed to calculate vibrational energy: $(e)")
    end

    return rho_evib[]
end

"""
$(SIGNATURES)

Calculate electron-electronic energy from state variables.

# Arguments
- `teex::Float64`: Electron-electronic temperature (K)
- `tvib::Float64`: Vibrational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)

# Returns
- `Float64`: Electron-electronic energy density

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
- `ArgumentError`: If teex ≤ 0, tvib ≤ 0, or arrays are invalid
"""
function calculate_electron_electronic_energy_wrapper(teex::Float64,
        tvib::Float64, rho_sp::Vector{Float64})
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    if teex <= 0.0
        throw(ArgumentError("Electron-electronic temperature must be positive, got: $teex"))
    end

    if tvib <= 0.0
        throw(ArgumentError("Vibrational temperature must be positive, got: $tvib"))
    end

    if length(rho_sp) == 0
        throw(ArgumentError("Species density array cannot be empty"))
    end

    rho_eeex = Ref{Float64}(0.0)

    # Call Fortran subroutine
    try
        ccall((:calculate_electron_electronic_energy, get_mtcr_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_eeex (output)
                Ref{Float64},                                    # teex
                Ref{Float64},                                    # tvib
                Ptr{Float64}),                                   # rho_sp
            rho_eeex,
            teex,
            tvib,
            rho_sp)
    catch e
        error("Failed to calculate electron-electronic energy: $(e)")
    end

    return rho_eeex[]
end

"""
$(SIGNATURES)

Set electronic state densities to Boltzmann distribution.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `tex::Float64`: Electronic temperature (K)
- `trot::Float64`: Rotational temperature (K)
- `tvib::Float64`: Vibrational temperature (K)

# Returns
- `Matrix{Float64}`: Electronic state densities in Boltzmann distribution

# Throws
- `ErrorException`: If MTCR library is not loaded or not initialized
- `ArgumentError`: If temperatures ≤ 0 or arrays are invalid
"""
function set_electronic_boltzmann_wrapper(rho_sp::Vector{Float64},
        tex::Float64, trot::Float64, tvib::Float64)
    if !is_mtcr_loaded()
        error("MTCR library not loaded. Call load_mtcr_library!(path) first.")
    end

    if !MTCR_INITIALIZED[]
        error("MTCR not initialized. Call initialize_api_wrapper() first.")
    end

    if tex <= 0.0 || trot <= 0.0 || tvib <= 0.0
        throw(ArgumentError("All temperatures must be positive. Got: tex=$tex, trot=$trot, tvib=$tvib"))
    end

    if length(rho_sp) == 0
        throw(ArgumentError("Species density array cannot be empty"))
    end

    # Get dimensions for electronic state array
    max_species = get_max_number_of_species_wrapper()
    # Get max atomic electronic states for rho_ex array
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()

    # Prepare output array
    rho_ex = zeros(Float64, max_atomic_electronic_states, max_species)

    # Call Fortran subroutine
    try
        ccall((:set_electronic_boltzmann, get_mtcr_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_ex (output)
                Ptr{Float64},                                    # rho_sp
                Ref{Float64},                                    # tex
                Ref{Float64},                                    # trot
                Ref{Float64}),                                  # tvib
            rho_ex,
            rho_sp,
            tex,
            trot,
            tvib)
    catch e
        error("Failed to set electronic Boltzmann distribution: $(e)")
    end

    return rho_ex
end
