
@testset "State Vector Layout" begin
    # Ensure Fortran-side library is loaded and initialized for dimension queries
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    # Use canonical nitrogen config for consistent dimensions
    config = mtcr.nitrogen_10ev_config()
    dims = mtcr.get_state_dimensions(config)

    # Basic dimension sanity
    @test dims.n_species == length(config.species)
    @test prod(dims.rho_ex_size) >= 0
    @test prod(dims.rho_vx_size) >= 0

    # Build consistent arrays and verify pack/unpack round-trip
    rho_sp = collect(1.0:(1.0 + dims.n_species - 1))
    rho_etot = 5.0
    rho_ex = fill(0.1, dims.rho_ex_size)
    rho_vx = fill(0.01, dims.rho_vx_size)

    u = mtcr.pack_state_vector(rho_sp, rho_etot, dims;
        rho_ex = rho_ex, rho_vx = rho_vx,
        rho_u = 0.0, rho_v = 0.0, rho_w = 0.0,
        rho_erot = 0.0, rho_eeex = 0.0, rho_evib = 0.0)

    expected_len = dims.n_species + 1 + prod(dims.rho_ex_size) + prod(dims.rho_vx_size) + 6
    @test length(u) == expected_len

    state = mtcr.unpack_state_vector(u, dims)
    @test state.rho_sp ≈ rho_sp
    @test state.rho_etot ≈ rho_etot
    @test state.rho_ex ≈ rho_ex
    @test state.rho_vx ≈ rho_vx
    @test state.rho_u == 0.0 && state.rho_v == 0.0 && state.rho_w == 0.0
    @test state.rho_erot == 0.0 && state.rho_eeex == 0.0 && state.rho_evib == 0.0

    # Derivative packing with and without optional arrays
    du = zeros(length(u))
    derivs = (
        drho_sp = ones(dims.n_species) * 0.1,
        drho_etot = 0.2,
        drho_ex = nothing,
        drho_vx = nothing,
        drho_erot = nothing,
        drho_eeex = nothing,
        drho_evib = nothing
    )
    @test_nowarn mtcr.pack_derivative_vector!(du, derivs, dims)
    @test du[1:(dims.n_species)] ≈ derivs.drho_sp
    @test du[dims.n_species + 1] ≈ derivs.drho_etot

    derivs_full = (
        drho_sp = ones(dims.n_species) * 0.1,
        drho_etot = 0.2,
        drho_ex = fill(0.01, dims.rho_ex_size),
        drho_vx = fill(0.001, dims.rho_vx_size),
        drho_erot = 0.05,
        drho_eeex = 0.02,
        drho_evib = 0.03
    )
    @test_nowarn mtcr.pack_derivative_vector!(du, derivs_full, dims)
    @test du[end - 2] ≈ 0.05
    @test du[end - 1] ≈ 0.02
    @test du[end] ≈ 0.03
end

@testset "RHS and Temperatures (initialized Fortran)" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    config = mtcr.nitrogen_10ev_config()
    state = mtcr.config_to_initial_state(config)
    dims = mtcr.get_state_dimensions(config)

    # Pack current state
    u = mtcr.pack_state_vector(state.rho_sp, state.rho_etot, dims;
        rho_ex = state.rho_ex,
        rho_eeex = state.rho_eeex,
        rho_evib = state.rho_evib)
    du = zeros(length(u))
    p = (
        dimensions = dims,
        config = config,
        rho_etot0 = state.rho_etot,
        molecular_weights = mtcr.get_molecular_weights(config.species)
    )

    # Compute RHS; expect finite derivatives
    @test_nowarn mtcr.mtcr_ode_system!(du, u, p, 0.0)
    @test all(isfinite, du)
    @test length(du) == length(u)

    # Temperature calculation at current state (charge-balanced electrons optional)
    temps = mtcr.calculate_temperatures_wrapper(state.rho_sp, state.rho_etot;
        rho_ex = state.rho_ex, rho_eeex = state.rho_eeex, rho_evib = state.rho_evib)
    @test temps.tt > 0 && temps.teex > 0 && temps.tvib > 0
end

@testset "Integrate 0D (solver pipeline)" begin
    # Initialize using the config-driven input to ensure the selected
    # database and options are honored (rather than a stale case file).

    config = mtcr.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    # Shorten integration to keep tests fast
    config = mtcr.MTCRConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = mtcr.TimeIntegrationConfig(5e-12, 1e-6, 1e-8, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        library_path = config.library_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_mtcr = false,
        print_source_terms = false
    )

    # Initialize the Fortran API using a temporary case generated from this config
    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = mtcr.config_to_initial_state(config)
    results = @time mtcr.integrate_0d_system(config, initial_state)
    @test results.time[end] > results.time[1]
    @test size(results.species_densities, 1) == length(config.species)
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test all(isfinite, results.temperatures.tv)
end

@testset "End-to-end Example (0D Adiabatic Nitrogen 10eV)" begin
    # Run the high-level example wrapper and verify structure and success
    results = @time mtcr.nitrogen_10ev_example()
    @test results.success == true
    @test length(results.time) >= 2
    @test size(results.species_densities, 1) == 5
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test all(isfinite, results.temperatures.tv)
    @test mtcr.validate_results(results)

    # Approximate final temperature values (update as needed)
    @test results.temperatures.tt[end]≈750.13 rtol=0.03
    @test results.temperatures.tv[end]≈758.49 rtol=0.03
    @test results.temperatures.te[end]≈2300.30 rtol=0.03

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.341e-14 rtol=0.05 # N
    @test results.species_densities[2, end]≈4.650e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈4.889e-19 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈4.760e-14 rtol=0.05 # N₂⁺
    @test results.species_densities[5, end]≈9.322e-19 rtol=0.05 # E⁻
end

# @testset "End-to-end Example (0D Isothermal Nitrogen 10eV)" begin
#     # Run the high-level example wrapper and verify structure and success
#     # results = @time mtcr.nitrogen_10ev_example(
#     #     test_case_path, isothermal = true)

#     results = @time mtcr.nitrogen_10ev_example(isothermal = true)
#     @test results.success == true
#     @test length(results.time) >= 2
#     @test size(results.species_densities, 1) == 5
#     @test all(isfinite, results.temperatures.tt)
#     @test all(isfinite, results.temperatures.te)
#     @test all(isfinite, results.temperatures.tv)
#     @test mtcr.validate_results(results)

#     # Approximate final temperature values (update as needed)
#     @test results.temperatures.tt[end]≈750.37 rtol=0.03
#     @test results.temperatures.tv[end]≈760.49 rtol=0.03
#     @test results.temperatures.te[end]≈2516.69 rtol=0.03

#     # Approximate final species densities in CGS (update as needed)
#     @test results.species_densities[1, end]≈8.921e-15 rtol=0.03 # N
#     @test results.species_densities[2, end]≈4.650e-10 rtol=0.03 # N₂
#     @test results.species_densities[3, end]≈2.635e-19 rtol=0.03 # N⁺
#     @test results.species_densities[4, end]≈4.946e-14 rtol=0.03 # N₂⁺
#     @test results.species_densities[5, end]≈9.685e-19 rtol=0.03 # E⁻
# end
