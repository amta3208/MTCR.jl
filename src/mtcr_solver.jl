"""
# MTCR Solver Module

This module provides the high-level interface for running MTCR simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.
"""

# Debug: RHS call counter
const _MTCR_ODE_DBG_CALLS = Ref(0)

"""
$(SIGNATURES)

Initialize the MTCR system.

This function must be called before any MTCR calculations can be performed.
It sets up the Fortran library, initializes internal data structures, and
prepares the system for simulation. The MTCR shared library is obtained from
the `MTCR_LIB_PATH` environment variable when it is not already loaded.

# Arguments
- `config::MTCRConfig`: Configuration for initialization
- `case_path::String`: Case directory path (optional, defaults to config.case_path)

# Returns
- `true` if initialization successful

# Throws
- `ErrorException` if initialization fails
"""
function initialize_mtcr(config::MTCRConfig, case_path::String = config.case_path)
    try
        # Ensure the shared library is loaded
        if !is_mtcr_loaded()
            load_mtcr_library!()
            @debug "MTCR library loaded successfully via MTCR_LIB_PATH"
        end

        # If the Fortran API is already initialized, finalize it so we can
        # reinitialize using the updated configuration and input files.
        try
            if is_api_initialized_wrapper()
                @debug "Finalizing existing MTCR API before re-initialization"
                finalize_api_wrapper()
            end
        catch e
            @warn "Unable to query/finalize existing MTCR API state" exception=e
        end

        # Generate input files for MTCR based on the provided configuration
        generate_input_files(config, case_path)
        @debug "MTCR input files generated" case_path=case_path

        # Initialize the API - get dimensions from Fortran
        result = initialize_api_wrapper(case_path = case_path)
        num_species = result.num_species
        num_dimensions = result.num_dimensions

        # Hard consistency check: configured species count must match MTCR setup
        if num_species != length(config.species)
            error("Configured species count ($(length(config.species))) does not match MTCR setup ($num_species). " *
                  "Ensure configuration and generated input match the MTCR database.")
        end

        @debug "MTCR initialized successfully" num_species=num_species num_dimensions=num_dimensions
        # Fetch and log runtime flags from MTCR for verification
        try
            flags = get_runtime_flags()
            @debug "MTCR runtime flags" ev_relax_set=flags.ev_relax_set vib_noneq=flags.vib_noneq eex_noneq=flags.eex_noneq rot_noneq=flags.rot_noneq bfe=flags.consider_elec_bfe bbh=flags.consider_elec_bbh bfh=flags.consider_elec_bfh bbe=flags.consider_elec_bbe
        catch e
            @warn "Unable to read MTCR runtime flags (rebuild library to enable)" exception=e
        end
        return true

    catch e
        @error "Failed to initialize MTCR" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Finalize the MTCR system and clean up resources.

This function should be called when MTCR is no longer needed to properly
clean up memory and resources.
"""
function finalize_mtcr()
    try
        # Call the wrapper finalization
        finalize_api_wrapper()

        # Close the library
        close_mtcr_library()

        @info "MTCR finalized successfully"
    catch e
        @error "Error during MTCR finalization" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Check if MTCR is properly initialized.

# Returns
- `true` if MTCR is initialized and ready for use
"""
function is_mtcr_initialized()
    try
        return is_mtcr_loaded() && is_api_initialized_wrapper()
    catch
        return false
    end
end

"""
$(SIGNATURES)

Convert MTCRConfig to initial state vectors for MTCR.

# Arguments
- `config::MTCRConfig`: Configuration object

# Returns
- Named tuple with initial state vectors in CGS units (as required by MTCR)
"""
function config_to_initial_state(config::MTCRConfig)
    # Get molecular weights
    molecular_weights = get_molecular_weights(config.species)

    # Ensure we're working in CGS - convert if needed
    config_cgs = config.unit_system == :CGS ? config : convert_config_units(config, :CGS)

    # Convert mole fractions to mass densities (CGS units)
    mass_densities_cgs = mole_fractions_to_mass_densities(
        config_cgs.mole_fractions, molecular_weights, config_cgs.total_number_density
    )

    # Calculate electronic state populations using MTCR's Boltzmann distribution
    # Use appropriate temperatures for each mode
    initial_electronic_states = set_electronic_boltzmann_wrapper(
        mass_densities_cgs,
        # Initialize electronic-state-resolved populations using the
        # electron-electronic temperature (TEE). Species that are not
        # electronically resolved do not contribute here (their mex = 0),
        # so using TEE only affects resolved species as intended.
        config_cgs.temperatures.Tee,  # Electronic-state populations use TEE
        config_cgs.temperatures.Tt,  # Rotational temperature (use translational as proxy)
        config_cgs.temperatures.Tv   # Vibrational temperature
    )

    # Calculate electron-electronic energy using MTCR's method.
    # IMPORTANT: Species without electronically resolved states contribute to
    # the electron-electronic mode; their energy must be initialized using the
    # electron temperature (TE), not TEE. MTCR internally excludes species with
    # resolved electronic states from this mode, so passing TE here is correct.
    initial_electron_electronic_energy = calculate_electron_electronic_energy_wrapper(
        config_cgs.temperatures.Te, config_cgs.temperatures.Tv, mass_densities_cgs
    )

    # Initialize vibrational state populations using MTCR's Boltzmann
    # distribution (for potential diagnostics), but do NOT use this to infer
    # STS. The active STS setting comes from the database/setup, which the
    # wrapper cannot reliably query here. To avoid a mismatch with MTCR,
    # always treat vibrational energy in mode form and do not include rho_vx
    # in Etot.
    initial_vibrational_states = set_vibrational_boltzmann_wrapper(
        initial_electronic_states,
        config_cgs.temperatures.Te,
        config_cgs.temperatures.Tt,
        config_cgs.temperatures.Tv
    )

    # Mode-level vibrational energy at Tv. Provide per-species electronic
    # temperature proxy as TE uniformly; species with resolved electronic
    # states are already accounted through rho_ex above.
    initial_vibrational_energy = calculate_vibrational_energy_wrapper(
        config_cgs.temperatures.Tv, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        tex = fill(config_cgs.temperatures.Te, length(mass_densities_cgs))
    )

    # Total energy: do not include rho_vx to avoid double counting when MTCR
    # is not in vibrational STS mode.
    initial_total_energy = calculate_total_energy_wrapper(
        config_cgs.temperatures.Tt, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        u = 0.0, v = 0.0, w = 0.0,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy
    )

    initial_remainder_energy = initial_total_energy - initial_electron_electronic_energy

    return (
        rho_sp = mass_densities_cgs,
        rho_etot = initial_total_energy,
        rho_rem = initial_remainder_energy,
        rho_ex = initial_electronic_states,
        rho_vx = nothing,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy,
        number_density = config_cgs.total_number_density,
        molecular_weights = molecular_weights,
        teex_const = config_cgs.temperatures.Te
    )
end

"""
$(SIGNATURES)

Calculate dimensions for the ODE state vector components.

# Arguments
- `config::MTCRConfig`: Configuration object

# Returns
- Named tuple with dimensions for each component
"""
function get_state_dimensions(config::MTCRConfig)
    if !is_mtcr_loaded()
        error("MTCR library must be loaded to get state dimensions. Set MTCR_LIB_PATH or call load_mtcr_library!(path) first.")
    end

    n_species = length(config.species)

    # Get actual dimensions from MTCR API
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()

    return (
        n_species = n_species,
        n_atomic_electronic_states = max_atomic_electronic_states,
        n_molecular_electronic_states = max_molecular_electronic_states,
        n_vibrational_levels = max_vibrational_quantum_number,
        rho_ex_size = (max_atomic_electronic_states, n_species),
        rho_vx_size = (
            max_vibrational_quantum_number + 1, max_molecular_electronic_states, n_species)  # +1 for 0:mnv indexing
    )
end

"""
$(SIGNATURES)

Compute the expected state vector length for a given dimensions tuple.

# Arguments
- `dimensions`: Dimensions structure from `get_state_dimensions()`

# Returns
- `Int`: Expected length of packed state/derivative vectors
"""
function expected_state_length(dimensions)
    n_species = dimensions.n_species
    nexc = prod(dimensions.rho_ex_size)
    nvib = prod(dimensions.rho_vx_size)
    # Layout: species + etot + ex(flat) + vx(flat) + (rho_u,rho_v,rho_w) + (rho_erot,rho_eeex,rho_evib)
    return n_species + 1 + nexc + nvib + 3 + 3
end

"""
$(SIGNATURES)

Pack state components into a single ODE state vector.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities
- `rho_etot::Float64`: Total energy density
- `dimensions`: Dimensions structure from get_state_dimensions()
- Optional components (all default to zero if not provided)

# Returns
- `Vector{Float64}`: Packed state vector
"""
function pack_state_vector(rho_sp::Vector{Float64}, rho_etot::Float64, dimensions;
        rho_ex = nothing, rho_vx = nothing, rho_u = nothing, rho_v = nothing,
        rho_w = nothing, rho_erot = nothing, rho_eeex = nothing, rho_evib = nothing)
    n_species = dimensions.n_species

    # Validate input array size
    if length(rho_sp) != n_species
        throw(ArgumentError("Species density array length ($(length(rho_sp))) must match number of species ($(n_species))"))
    end

    # Preallocate output for determinism and speed
    u = Vector{Float64}(undef, expected_state_length(dimensions))
    idx = 1

    # Species densities
    @inbounds begin
        u[idx:(idx + n_species - 1)] .= rho_sp
    end
    idx += n_species

    # Total energy density
    u[idx] = rho_etot
    idx += 1

    # Electronic state densities (flattened)
    ex_target = dimensions.rho_ex_size
    if rho_ex !== nothing
        if size(rho_ex) != ex_target
            ex_fixed = zeros(Float64, ex_target)
            m1 = min(size(rho_ex, 1), ex_target[1])
            m2 = min(size(rho_ex, 2), ex_target[2])
            @inbounds ex_fixed[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
            @inbounds u[idx:(idx + prod(ex_target) - 1)] .= vec(ex_fixed)
        else
            @inbounds u[idx:(idx + prod(ex_target) - 1)] .= vec(rho_ex)
        end
    else
        @inbounds fill!(view(u, idx:(idx + prod(ex_target) - 1)), 0.0)
    end
    idx += prod(ex_target)

    # Vibrational state densities (flattened)
    vx_target = dimensions.rho_vx_size
    if rho_vx !== nothing
        if size(rho_vx) != vx_target
            vx_fixed = zeros(Float64, vx_target)
            m1 = min(size(rho_vx, 1), vx_target[1])
            m2 = min(size(rho_vx, 2), vx_target[2])
            m3 = min(size(rho_vx, 3), vx_target[3])
            @inbounds vx_fixed[1:m1, 1:m2, 1:m3] .= rho_vx[1:m1, 1:m2, 1:m3]
            @inbounds u[idx:(idx + prod(vx_target) - 1)] .= vec(vx_fixed)
        else
            @inbounds u[idx:(idx + prod(vx_target) - 1)] .= vec(rho_vx)
        end
    else
        @inbounds fill!(view(u, idx:(idx + prod(vx_target) - 1)), 0.0)
    end
    idx += prod(vx_target)

    # Velocity components
    u[idx] = (rho_u === nothing ? 0.0 : rho_u)
    idx += 1
    u[idx] = (rho_v === nothing ? 0.0 : rho_v)
    idx += 1
    u[idx] = (rho_w === nothing ? 0.0 : rho_w)
    idx += 1

    # Energy mode components
    u[idx] = (rho_erot === nothing ? 0.0 : rho_erot)
    idx += 1
    u[idx] = (rho_eeex === nothing ? 0.0 : rho_eeex)
    idx += 1
    u[idx] = (rho_evib === nothing ? 0.0 : rho_evib)

    # Final sanity check
    @assert idx==length(u) "Internal error: packed state length mismatch"

    return u
end

"""
$(SIGNATURES)

Unpack ODE state vector into individual components.

# Arguments
- `u::Vector{Float64}`: Packed state vector
- `dimensions`: Dimensions structure from get_state_dimensions()

# Returns
- Named tuple with unpacked components
"""
function unpack_state_vector(u::Vector{Float64}, dimensions)
    n_species = dimensions.n_species
    rho_ex_size = dimensions.rho_ex_size
    rho_vx_size = dimensions.rho_vx_size

    idx = 1

    # Extract species densities (return a concrete Vector for wrapper calls)
    rho_sp = Vector{Float64}(@view u[idx:(idx + n_species - 1)])
    idx += n_species

    # Extract total energy density
    rho_etot = u[idx]
    idx += 1

    # Extract electronic state densities
    rho_ex_flat = @view u[idx:(idx + prod(rho_ex_size) - 1)]
    rho_ex = reshape(Vector{Float64}(rho_ex_flat), rho_ex_size)
    idx += prod(rho_ex_size)

    # Extract vibrational state densities
    rho_vx_flat = @view u[idx:(idx + prod(rho_vx_size) - 1)]
    rho_vx = reshape(Vector{Float64}(rho_vx_flat), rho_vx_size)
    idx += prod(rho_vx_size)

    # Extract tail components (rho_u, rho_v, rho_w, rho_erot, rho_eeex, rho_evib)
    rho_u = u[idx]
    idx += 1
    rho_v = u[idx]
    idx += 1
    rho_w = u[idx]
    idx += 1
    rho_erot = u[idx]
    idx += 1
    rho_eeex = u[idx]
    idx += 1
    rho_evib = u[idx]

    return (
        rho_sp = rho_sp,
        rho_etot = rho_etot,
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_u = rho_u,
        rho_v = rho_v,
        rho_w = rho_w,
        rho_erot = rho_erot,
        rho_eeex = rho_eeex,
        rho_evib = rho_evib
    )
end

"""
$(SIGNATURES)

Pack derivative components into ODE derivative vector.

# Arguments
- `du::Vector{Float64}`: Output derivative vector (modified in-place)
- `derivatives`: Derivatives from calculate_sources_wrapper()
- `dimensions`: Dimensions structure

# Returns
- Nothing (modifies du in-place)
"""
function pack_derivative_vector!(du::Vector{Float64}, derivatives, dimensions)
    n_species = dimensions.n_species
    rho_ex_size = dimensions.rho_ex_size
    rho_vx_size = dimensions.rho_vx_size

    idx = 1

    # Pack species density derivatives
    du[idx:(idx + n_species - 1)] .= derivatives.drho_sp
    idx += n_species

    # Pack total energy derivative
    du[idx] = derivatives.drho_etot
    idx += 1

    # Pack electronic state derivatives
    if derivatives.drho_ex !== nothing
        if size(derivatives.drho_ex) != rho_ex_size
            throw(DimensionMismatch("drho_ex size $(size(derivatives.drho_ex)) does not match expected $(rho_ex_size)"))
        end
        du[idx:(idx + prod(rho_ex_size) - 1)] .= vec(derivatives.drho_ex)
    else
        du[idx:(idx + prod(rho_ex_size) - 1)] .= 0.0
    end
    idx += prod(rho_ex_size)

    # Pack vibrational state derivatives
    if derivatives.drho_vx !== nothing
        if size(derivatives.drho_vx) != rho_vx_size
            throw(DimensionMismatch("drho_vx size $(size(derivatives.drho_vx)) does not match expected $(rho_vx_size)"))
        end
        du[idx:(idx + prod(rho_vx_size) - 1)] .= vec(derivatives.drho_vx)
    else
        du[idx:(idx + prod(rho_vx_size) - 1)] .= 0.0
    end
    idx += prod(rho_vx_size)

    # Pack velocity derivatives (zero for 0D)
    du[idx] = 0.0  # drho_u
    idx += 1
    du[idx] = 0.0  # drho_v
    idx += 1
    du[idx] = 0.0  # drho_w
    idx += 1

    # Pack energy derivatives
    du[idx] = derivatives.drho_erot === nothing ? 0.0 : derivatives.drho_erot
    idx += 1
    du[idx] = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
    idx += 1
    du[idx] = derivatives.drho_evib === nothing ? 0.0 : derivatives.drho_evib

    # Sanity check: ensure we've filled the entire vector
    @assert idx==length(du) "Internal error: derivative vector length mismatch"

    return nothing
end

"""
Reconstruct energy components based on configuration flags.

When the isothermal electron-electronic mode is enabled, the energy slot stored
in the state vector holds the remainder energy (`rho_rem`). This helper
recomputes the electron-electronic energy from the prescribed `Tee`, rebuilds
the total energy, and returns all relevant pieces for downstream use. When the
mode is disabled, it simply returns the existing state components.
"""
function reconstruct_energy_components(state, config;
        teex_const::Float64 = config.temperatures.Te,
        teex_vec::Union{Nothing, AbstractVector{<:Real}} = nothing)
    is_isothermal = config.physics.is_isothermal_teex

    if !is_isothermal
        rho_etot = state.rho_etot
        rho_eeex = state.rho_eeex
        rho_rem = rho_etot - rho_eeex
        return (rho_etot = rho_etot,
            rho_eeex = rho_eeex,
            rho_rem = rho_rem,
            tvib = nothing)
    end

    rho_rem = state.rho_etot
    rho_sp = state.rho_sp

    teex_kw = teex_vec === nothing ? nothing : teex_vec

    tvib = calculate_vibrational_temperature_wrapper(
        state.rho_evib, rho_sp;
        rho_ex = state.rho_ex,
        tex = teex_kw)

    rho_eeex = calculate_electron_electronic_energy_wrapper(teex_const, tvib, rho_sp)
    rho_etot = rho_rem + rho_eeex

    return (rho_etot = rho_etot,
        rho_eeex = rho_eeex,
        rho_rem = rho_rem,
        tvib = tvib)
end

"""
$(SIGNATURES)

ODE system function for MTCR integration.

This function defines the ODE system du/dt = f(u, p, t) where:
- u is the state vector containing all MTCR variables
- p contains parameters (dimensions, etc.)
- t is time

# Arguments
- `du::Vector{Float64}`: Output derivative vector
- `u::Vector{Float64}`: Input state vector
- `p`: Parameters structure
- `t::Float64`: Current time

# Returns
- Nothing (modifies du in-place)
"""
function mtcr_ode_system!(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    # Guard against invalid inputs to the Fortran layer
    if any(!isfinite, u)
        fill!(du, 0.0)
        return nothing
    end
    nsp = p.dimensions.n_species
    if any(@view(u[1:nsp]) .< 0.0)
        fill!(du, 0.0)
        return nothing
    end
    if u[nsp + 1] < 0.0
        fill!(du, 0.0)
        return nothing
    end
    # Unpack state vector
    state = unpack_state_vector(u, p.dimensions)

    # Calculate source terms using MTCR
    try
        config = hasproperty(p, :config) ? p.config : nothing
        teex_vec = hasproperty(p, :teex_const_vec) ? p.teex_const_vec : nothing
        teex_const = hasproperty(p, :teex_const) ? p.teex_const :
                     (config === nothing ? 0.0 : config.temperatures.Te)
        is_isothermal = config !== nothing && config.physics.is_isothermal_teex

        energy = config === nothing ?
                 (rho_etot = state.rho_etot,
            rho_eeex = state.rho_eeex,
            rho_rem = state.rho_etot - state.rho_eeex,
            tvib = nothing) :
                 reconstruct_energy_components(state, config;
            teex_const = teex_const,
            teex_vec = teex_vec)

        rho_etot_effective = energy.rho_etot
        rho_eeex_effective = energy.rho_eeex

        # Always pass the current total energy to the Fortran RHS so that
        # temperatures and rates are computed from the instantaneous state.
        # Mirror MTCR handling of total energy:
        # - If radiation is OFF, hold total energy constant across RHS calls
        # - If radiation is ON, allow MTCR to evolve total energy
        use_const_etot = !is_isothermal && hasproperty(p, :config) &&
                         hasproperty(p.config, :processes) &&
                         getfield(p.config.processes, :consider_rad) == 0
        rho_etot_to_use = use_const_etot && hasproperty(p, :rho_etot0) ? p.rho_etot0 :
                          rho_etot_effective

        derivatives = calculate_sources_wrapper(
            state.rho_sp, rho_etot_to_use;
            rho_ex = state.rho_ex,
            rho_vx = (sum(state.rho_vx) > 0.0) ? state.rho_vx : nothing,
            rho_u = state.rho_u,
            rho_v = state.rho_v,
            rho_w = state.rho_w,
            rho_erot = state.rho_erot,
            rho_eeex = rho_eeex_effective,
            rho_evib = state.rho_evib
        )

        # Optional solver-side debug mirror of wrapper outputs
        if get(ENV, "MTCR_SOLVER_DEBUG", "0") == "1"
            @info "SOLVER_DERIV" drho_etot=derivatives.drho_etot drho_erot=derivatives.drho_erot drho_eeex=derivatives.drho_eeex drho_evib=derivatives.drho_evib
        end

        # If solving electrons via charge balance, adjust electron derivative
        if hasproperty(p, :config) &&
           p.config.physics.get_electron_density_by_charge_balance
            names = p.config.species
            weights = p.molecular_weights
            # simple charge inference from species names
            charges = map(nm -> count(==('+'), nm) - count(==('-'), nm), names)
            elec_idx = findfirst(==("E-"), names)
            elec_idx = elec_idx === nothing ? findfirst(==(-1), charges) : elec_idx
            if elec_idx !== nothing
                spwt_e = weights[elec_idx]
                s = 0.0
                @inbounds for i in eachindex(derivatives.drho_sp)
                    if i == elec_idx
                        continue
                    end
                    zi = charges[i]
                    if zi != 0
                        s += (zi / weights[i]) * derivatives.drho_sp[i]
                    end
                end
                new_drho_sp = copy(derivatives.drho_sp)
                new_drho_sp[elec_idx] = spwt_e * s
                derivatives = (
                    drho_sp = new_drho_sp,
                    drho_etot = derivatives.drho_etot,
                    drho_ex = derivatives.drho_ex,
                    drho_vx = derivatives.drho_vx,
                    drho_erot = derivatives.drho_erot,
                    drho_eeex = derivatives.drho_eeex,
                    drho_evib = derivatives.drho_evib
                )
            end
        end

        if is_isothermal
            drho_eeex_val = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
            drho_etot_raw = derivatives.drho_etot
            drho_rem = drho_etot_raw - drho_eeex_val
            derivatives = merge(derivatives, (drho_etot = drho_rem, drho_eeex = 0.0))
        end

        # Debug: print derivative norms on first few RHS calls (debug-level only)
        _MTCR_ODE_DBG_CALLS[] += 1
        if _MTCR_ODE_DBG_CALLS[] <= 5
            drsp_norm = maximum(abs, derivatives.drho_sp)
            dr_etot = derivatives.drho_etot
            dr_eeex = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
            dr_evib = derivatives.drho_evib === nothing ? 0.0 : derivatives.drho_evib
            @debug "MTCR RHS debug" t=t drho_sp_max=drsp_norm drho_etot=dr_etot drho_eeex=dr_eeex drho_evib=dr_evib
        end

        # Pack derivatives into output vector
        pack_derivative_vector!(du, derivatives, p.dimensions)
        if use_const_etot
            n_species = p.dimensions.n_species
            du[n_species + 1] = 0.0
        end

    catch e
        @error "Error in MTCR ODE system" exception=e time=t
        # Fill with zeros to prevent integration failure
        fill!(du, 0.0)
    end

    return nothing
end

"""
$(SIGNATURES)

Integrate the 0D system over time using DifferentialEquations.jl.

# Arguments
- `config::MTCRConfig`: Configuration object
- `initial_state`: Initial state vectors in CGS units

# Returns
- `MTCRResults`: Simulation results (converted back to SI units)
"""
function integrate_0d_system(config::MTCRConfig, initial_state)
    # Time parameters
    dt = config.time_params.dt
    dtm = config.time_params.dtm
    tlim = config.time_params.tlim
    nstep = config.time_params.nstep

    @info "Setting up ODE integration" tlim=tlim

    # Get state dimensions
    dimensions = get_state_dimensions(config)
    n_species = dimensions.n_species
    is_isothermal = config.physics.is_isothermal_teex

    # Create initial ODE state vector (all in CGS units)
    # Include available electronic states and energy components so the ODE RHS
    # receives consistent nonequilibrium energies at t=0.
    energy_scalar0 = is_isothermal ? initial_state.rho_rem : initial_state.rho_etot
    rho_eeex0 = is_isothermal ? 0.0 : initial_state.rho_eeex
    u0 = pack_state_vector(
        initial_state.rho_sp, energy_scalar0, dimensions;
        rho_ex = initial_state.rho_ex,
        # rho_vx = initial_state.rho_vx,
        # rho_eeex = initial_state.rho_eeex,
        rho_eeex = rho_eeex0,
        rho_evib = initial_state.rho_evib
    )

    @info "Initial state vector created" length_u0=length(u0) n_species=n_species

    # Time span for integration
    tspan = (0.0, tlim)

    # Parameters for the ODE system. Keep total energy density constant
    # (aligned with Fortran API behavior) by passing it as a parameter.
    molecular_weights = get_molecular_weights(config.species)
    teex_const_vec = fill(config.temperatures.Te, n_species)
    p = (dimensions = dimensions, config = config, rho_etot0 = initial_state.rho_etot,
        molecular_weights = molecular_weights,
        teex_const = initial_state.teex_const,
        teex_const_vec = teex_const_vec)

    # Create ODE problem (we add domain-safety via a callback below)
    prob = ODEProblem(mtcr_ode_system!, u0, tspan, p)

    @info "ODE problem created, starting integration..."

    try
        sol = solve(prob;
            alg_hints = [:stiff],
            dt = dt,
            # dtmax = dtm,
            reltol = 1e-11,
            abstol = 1e-13,
            save_everystep = true,
            verbose = true
        )

        @info "ODE integration completed" success=sol.retcode

        # MTCR-style status printing for each solver step
        iter = 0
        # Use the actually-taken first step size when available
        first_dt = length(sol.t) >= 2 ? (sol.t[2] - sol.t[1]) : dt
        for (i, t) in enumerate(sol.t)
            st = unpack_state_vector(sol.u[i], dimensions)
            energy = reconstruct_energy_components(st, config;
                teex_const = config.temperatures.Te,
                teex_vec = teex_const_vec)
            # Charge-balanced electron density for diagnostics
            rho_sp_cb = copy(st.rho_sp)
            if config.physics.get_electron_density_by_charge_balance
                names = config.species
                charges = map(nm -> count(==('+'), nm) - count(==('-'), nm), names)
                elec_idx = findfirst(==("E-"), names)
                elec_idx = elec_idx === nothing ? findfirst(==(-1), charges) : elec_idx
                if elec_idx !== nothing
                    spwt_e = molecular_weights[elec_idx]
                    s = 0.0
                    @inbounds for k in eachindex(rho_sp_cb)
                        if k == elec_idx
                            continue
                        end
                        zk = charges[k]
                        if zk != 0
                            s += (zk / molecular_weights[k]) * rho_sp_cb[k]
                        end
                    end
                    rho_sp_cb[elec_idx] = spwt_e * s
                end
            end

            temps = calculate_temperatures_wrapper(rho_sp_cb, energy.rho_etot;
                rho_ex = st.rho_ex,
                rho_eeex = energy.rho_eeex, rho_evib = st.rho_evib)

            # Mole fractions
            denom = sum(rho_sp_cb ./ molecular_weights)
            x = (rho_sp_cb ./ molecular_weights) ./ denom

            # Mass fraction sum error
            ys = st.rho_sp ./ sum(st.rho_sp)
            yerr = abs(sum(ys) - 1.0)
            println(@sprintf(" Ytot,err   = % .3E", yerr))

            # Relative enthalpy change (%) at current Tt vs reference energy
            Ecomp = calculate_total_energy_wrapper(temps.tt, rho_sp_cb;
                rho_ex = st.rho_ex,
                rho_eeex = energy.rho_eeex, rho_evib = st.rho_evib)
            dEnth = 100.0 * (Ecomp - energy.rho_etot) / energy.rho_etot
            println(@sprintf(" dEnth (%%)  = % .5E", dEnth))

            t_us = t * 1e6
            dt_us = (i == 1 ? first_dt : (sol.t[i] - sol.t[i - 1])) * 1e6
            println(@sprintf(" iter       = %6d", iter))
            println(@sprintf(" time       = % .2E mu-s ", t_us))
            println(@sprintf(" dt         = % .2E mu-s", dt_us))
            println(@sprintf(" T(t,e,r,v) = % .3E % .3E % .3E % .3E K",
                temps.tt, temps.teex, temps.trot, temps.tvib))

            # Species mole fractions line
            xbuf = IOBuffer()
            print(xbuf, " X          =")
            for xi in x
                print(xbuf, @sprintf(" % .3E", xi))
            end
            println(String(take!(xbuf)))
            println()

            iter += 1
        end

        # Extract results at output times
        n_times = length(sol.t)
        time_points = collect(sol.t)

        # Pre-allocate arrays
        species_densities = zeros(n_species, n_times)
        temperatures_tt = Vector{Float64}(undef, n_times)
        temperatures_te = Vector{Float64}(undef, n_times)
        temperatures_tv = Vector{Float64}(undef, n_times)
        total_energies = Vector{Float64}(undef, n_times)

        # Extract species densities and calculate temperatures at each time point
        for i in 1:n_times
            state = unpack_state_vector(sol.u[i], dimensions)
            energy = reconstruct_energy_components(state, config;
                teex_const = config.temperatures.Te,
                teex_vec = teex_const_vec)

            # Store species densities
            species_densities[:, i] = state.rho_sp
            total_energies[i] = energy.rho_etot

            # Calculate temperatures for this time point using charge-balanced electrons when requested
            rho_sp_cb = copy(state.rho_sp)
            if config.physics.get_electron_density_by_charge_balance
                names = config.species
                charges = map(nm -> count(==('+'), nm) - count(==('-'), nm), names)
                elec_idx = findfirst(==("E-"), names)
                elec_idx = elec_idx === nothing ? findfirst(==(-1), charges) : elec_idx
                if elec_idx !== nothing
                    spwt_e = molecular_weights[elec_idx]
                    s = 0.0
                    @inbounds for k in eachindex(rho_sp_cb)
                        if k == elec_idx
                            continue
                        end
                        zk = charges[k]
                        if zk != 0
                            s += (zk / molecular_weights[k]) * rho_sp_cb[k]
                        end
                    end
                    rho_sp_cb[elec_idx] = spwt_e * s
                end
            end

            # Only pass rho_vx when it has nonzero content (STS active)
            rho_vx_arg = (sum(state.rho_vx) > 0) ? state.rho_vx : nothing
            try
                temps = calculate_temperatures_wrapper(
                    rho_sp_cb, energy.rho_etot; rho_ex = state.rho_ex, rho_vx = rho_vx_arg,
                    rho_eeex = energy.rho_eeex, rho_evib = state.rho_evib)

                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib
            catch e
                @warn "Temperature calculation failed at time $(time_points[i])" exception=e
                # Use previous values or reasonable defaults
                if i > 1
                    temperatures_tt[i] = temperatures_tt[i - 1]
                    temperatures_te[i] = temperatures_te[i - 1]
                    temperatures_tv[i] = temperatures_tv[i - 1]
                else
                    temperatures_tt[i] = config.temperatures.Tt
                    temperatures_te[i] = config.temperatures.Te
                    temperatures_tv[i] = config.temperatures.Tv
                end
            end
        end

        # Convert results back to SI units if needed
        if config.unit_system == :SI
            species_densities_si = zeros(size(species_densities))
            for i in axes(species_densities, 2)
                species_densities_si[:, i] = convert_density_cgs_to_si(species_densities[
                    :, i])
            end
            total_energies_si = [convert_energy_density_cgs_to_si(e)
                                 for e in total_energies]
        else
            species_densities_si = species_densities
            total_energies_si = total_energies
        end

        # Create results structure
        temperatures = (
            tt = temperatures_tt,
            te = temperatures_te,
            tv = temperatures_tv
        )

        # Robust success detection across SciMLBase versions
        rc = sol.retcode
        success = rc isa Symbol ? (rc in (:Success, :Terminated)) :
                  (occursin("Success", string(rc)) || occursin("Terminated", string(rc)))
        message = success ? "ODE integration completed successfully" :
                  "ODE integration terminated: $(rc)"

        @info "Results processing completed" final_time=time_points[end] success=success

        return MTCRResults(
            time_points,
            species_densities_si,
            temperatures,
            total_energies_si,
            nothing,  # source_terms - could be added later
            success,
            message
        )

    catch e
        @error "ODE integration failed" exception=e

        # Return minimal results with error information
        return MTCRResults(
            [0.0],
            reshape(initial_state.rho_sp, :, 1),
            (tt = [config.temperatures.Tt], te = [config.temperatures.Te],
                tv = [config.temperatures.Tv], tee = [config.temperatures.Te]),
            [initial_state.rho_etot],
            nothing,
            false,
            "ODE integration failed: $(string(e))"
        )
    end
end

"""
$(SIGNATURES)

Solve a 0D MTCR simulation.

This is the main high-level interface for running MTCR simulations.
It handles all the complexity of data conversion, Fortran interfacing,
and result processing.

# Arguments
- `config::MTCRConfig`: Configuration for the simulation

# Returns
- `MTCRResults`: Results of the simulation

# Throws
- `ErrorException` if MTCR not initialized or simulation fails
"""
function solve_mtcr_0d(config::MTCRConfig)
    if !is_mtcr_initialized()
        error("MTCR not initialized. Call initialize_mtcr(config) first.")
    end

    try
        @info "Starting MTCR 0D simulation" species=config.species

        # Convert configuration to initial conditions (SI to CGS)
        initial_state = config_to_initial_state(config)

        # Run the time integration
        results = integrate_0d_system(config, initial_state)

        @info "MTCR simulation completed successfully"
        return results

    catch e
        @error "MTCR simulation failed" exception=e
        return MTCRResults(
            Float64[], zeros(0, 0), (;), Float64[], nothing, false,
            "Simulation failed: $(string(e))"
        )
    end
end

"""
$(SIGNATURES)

Run the 0D Nitrogen Te=10eV example case.

This function provides a convenient way to run the reference test case
that matches the MTCR example in `/mtcr/examples/0D_Nitrogen_Te_10eV`.
Requires the `MTCR_LIB_PATH` environment variable to point to the MTCR shared library.

# Arguments
- `case_path::String`: Case directory path (optional, creates temp directory if not provided)

# Returns
- `MTCRResults`: Results of the simulation

# Example
```julia
results = nitrogen_10ev_example()
```
"""
function nitrogen_10ev_example(case_path::String = mktempdir();
        isothermal::Bool = false)
    # Create configuration for the example case
    config = nitrogen_10ev_config(; isothermal = isothermal)

    # Update config with case path
    config_with_path = MTCRConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = config.time_params,
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        library_path = config.library_path,
        case_path = case_path,
        unit_system = config.unit_system,
        validate_species_against_mtcr = config.validate_species_against_mtcr,
        print_source_terms = config.print_source_terms
    )

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config_with_path.species mole_fractions=config_with_path.mole_fractions
    @info "Temperatures" Tt=config_with_path.temperatures.Tt Te=config_with_path.temperatures.Te
    @info "Time parameters" dt=config_with_path.time_params.dt tlim=config_with_path.time_params.tlim
    @info "Case path" case_path=case_path

    try
        # Initialize MTCR with config
        initialize_mtcr(config_with_path, case_path)

        # Run simulation
        results = solve_mtcr_0d(config_with_path)

        if results.success
            @info "Example simulation completed successfully"
            @info "Final conditions" time=results.time[end]

            # Print final species densities
            for (i, species) in enumerate(config_with_path.species)
                final_density = results.species_densities[i, end]
                unit_str = config_with_path.unit_system == :SI ? "kg/m³" : "g/cm³"
                @info "Final density" species=species density=final_density unit=unit_str
            end

            @info "Final temperatures" Tt=results.temperatures.tt[end] Tv=results.temperatures.tv[end] Te=results.temperatures.te[end]
        else
            @error "Example simulation failed" message=results.message
        end

        return results

    finally
        # Clean up MTCR resources
        try
            finalize_mtcr()
        catch e
            @warn "Error during cleanup" exception=e
        end
    end
end

"""
$(SIGNATURES)

Validate simulation results for physical consistency.

# Arguments
- `results::MTCRResults`: Simulation results to validate

# Returns
- `true` if results pass validation, `false` otherwise
"""
function validate_results(results::MTCRResults)
    if !results.success
        @warn "Simulation was not successful"
        return false
    end

    # Check for negative densities
    if any(results.species_densities .< 0)
        @warn "Negative species densities found"
        return false
    end

    # Check for NaN or Inf values
    if any(isnan.(results.species_densities)) || any(isinf.(results.species_densities))
        @warn "NaN or Inf values found in species densities"
        return false
    end

    if any(isnan.(results.temperatures.tt)) || any(isinf.(results.temperatures.tt))
        @warn "NaN or Inf values found in temperatures"
        return false
    end

    # Check temperature ranges (should be positive and reasonable)
    if any(results.temperatures.tt .<= 0) || any(results.temperatures.tt .> 1e6)
        @warn "Unreasonable translational temperatures found"
        return false
    end

    if any(results.temperatures.te .<= 0) || any(results.temperatures.te .> 1e6)
        @warn "Unreasonable electron temperatures found"
        return false
    end

    @info "Results validation passed"
    return true
end

"""
$(SIGNATURES)

Save MTCR results to file.

# Arguments
- `results::MTCRResults`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::MTCRResults, filename::String)
    try
        # Prepare data for CSV output
        n_times = length(results.time)
        n_species = size(results.species_densities, 1)

        # Create header
        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        # Create data matrix
        data = zeros(n_times, length(header))
        data[:, 1] = results.time
        data[:, 2] = results.total_energy
        data[:, 3] = results.temperatures.tt
        data[:, 4] = results.temperatures.te
        data[:, 5] = results.temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = results.species_densities[i, :]
        end

        # Write to file
        open(filename, "w") do io
            # Write header
            println(io, join(header, ","))

            # Write data
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename=filename
        return true

    catch e
        @error "Failed to save results" filename=filename exception=e
        return false
    end
end
