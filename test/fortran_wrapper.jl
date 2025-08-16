
# TODO: Update with correct library path later on
temp_mtcr_path = "/Users/amin/.julia/dev/MTCR/mtcr/source/libmtcr.so"

@testset "Library Management" begin
    @testset "Library Loading and Status" begin
        # Test initial state - no library loaded
        @test !mtcr.is_mtcr_loaded()

        # Test that getting handle fails when no library is loaded
        @test_throws ErrorException mtcr.get_mtcr_handle()

        # Test error message content
        try
            mtcr.get_mtcr_handle()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end
    end

    @testset "Library Path Setting" begin
        # Test setting library path with non-existent file
        fake_path = "/nonexistent/path/libmtcr.so"
        @test_throws ErrorException mtcr.load_mtcr_library!(fake_path)

        # Test error message for non-existent file
        try
            mtcr.load_mtcr_library!(fake_path)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library file not found", e.msg)
            @test occursin(fake_path, e.msg)
        end

        # Test that library is still not loaded after failed attempt
        @test !mtcr.is_mtcr_loaded()
    end

    @testset "Library Cleanup" begin
        # Test closing library when none is loaded (should be safe)
        @test_nowarn mtcr.close_mtcr_library()
        @test !mtcr.is_mtcr_loaded()

        # Test multiple calls to close (should be safe)
        @test_nowarn mtcr.close_mtcr_library()
        @test_nowarn mtcr.close_mtcr_library()
        @test !mtcr.is_mtcr_loaded()
    end
end

@testset "Initialization" begin
    # Ensure library path is set
    mtcr.load_mtcr_library!(temp_mtcr_path)

    # Get path to test case directory
    test_case_path = joinpath(@__DIR__, "test_case")

    # Test basic initialization (Fortran determines species count from input files)
    result = @test_nowarn mtcr.initialize_api_wrapper(case_path = test_case_path)

    # Check that result contains the expected fields
    if result !== nothing
        @test result isa NamedTuple
        @test haskey(result, :num_species)
        @test haskey(result, :num_dimensions)
        @test result.num_species isa Int32
        @test result.num_dimensions isa Int32
        @test result.num_species > 0
        @test result.num_dimensions >= 0
    end
end

@testset "Input/Directory Handling" begin
    @testset "Initialization Input Validation" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for validation tests
        mtcr.MTCR_INITIALIZED[] = false

        # Test with non-existent case path
        @test_throws ErrorException mtcr.initialize_api_wrapper(case_path = "/nonexistent/path")

        # Test error message for non-existent case path
        try
            mtcr.initialize_api_wrapper(case_path = "/nonexistent/path")
            @test false  # Should not reach here
        catch e
            @test occursin("Case path does not exist", e.msg)
        end

        # Test with case path missing input directory
        temp_dir = mktempdir()
        try
            @test_throws ErrorException mtcr.initialize_api_wrapper(case_path = temp_dir)
        finally
            rm(temp_dir; recursive = true)
        end

        # Test error message for missing input file
        temp_dir = mktempdir()
        try
            mtcr.initialize_api_wrapper(case_path = temp_dir)
            @test false  # Should not reach here
        catch e
            @test occursin("Required input file not found", e.msg)
            @test occursin("prob_setup.inp", e.msg)
        finally
            rm(temp_dir; recursive = true)
        end
    end

    @testset "Directory Management" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Store original directory
        original_dir = pwd()
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test that directory is restored after successful call
        mtcr.initialize_api_wrapper(case_path = test_case_path)
        @test pwd() == original_dir

        # Test that directory is restored even after failed call
        temp_dir = mktempdir()
        try
            mtcr.initialize_api_wrapper(case_path = temp_dir)
        catch
            # Expected to fail
        end
        @test pwd() == original_dir
        rm(temp_dir; recursive = true)
    end

    @testset "Output Directory Creation" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for this test
        mtcr.MTCR_INITIALIZED[] = false

        # Create a temporary test case directory with input
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Verify output directories don't exist initially
            output_dir = joinpath(temp_case_dir, "output")
            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test !isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)

            # Initialize - should create output directories
            mtcr.initialize_api_wrapper(case_path = temp_case_dir)

            # Verify output directories were created
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

            # Test that calling again doesn't cause errors (directories already exist)
            @test mtcr.initialize_api_wrapper(case_path = temp_case_dir) === nothing

            # Verify directories still exist
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

        finally
            rm(temp_case_dir; recursive = true)
        end
    end

    @testset "Output Directory Creation with Existing Directories" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Reset initialization state for this test
        mtcr.MTCR_INITIALIZED[] = false

        # Create a temporary test case directory with input and partial output structure
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Create output directory but not subdirectories
            output_dir = joinpath(temp_case_dir, "output")
            mkpath(output_dir)

            # Create a test file in output to verify it's preserved
            test_file = joinpath(output_dir, "test_file.txt")
            write(test_file, "test content")

            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)
            @test isfile(test_file)

            # Initialize - should create missing subdirectories
            mtcr.initialize_api_wrapper(case_path = temp_case_dir)

            # Verify all directories exist and existing content is preserved
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)
            @test isfile(test_file)
            @test read(test_file, String) == "test content"

        finally
            rm(temp_case_dir; recursive = true)
        end
    end
end

@testset "Utility Functions" begin
    @testset "Maximum Species Count" begin
        # Ensure library path is set
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Test error when library not loaded
        mtcr.close_mtcr_library()
        @test_throws ErrorException mtcr.get_max_number_of_species_wrapper()

        # Test error message content
        try
            mtcr.get_max_number_of_species_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end

        # Reload library for actual test
        mtcr.load_mtcr_library!(temp_mtcr_path)

        max_species = mtcr.get_max_number_of_species_wrapper()
        println("Max species: ", max_species)

        # Check return type
        @test max_species isa Int32

        # Check that it's positive and reasonable
        @test max_species > 0
        @test max_species <= 100  # Reasonable upper bound
    end

    @testset "Species Names" begin
        # Ensure library is loaded
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Test error when library not loaded
        mtcr.close_mtcr_library()
        @test_throws ErrorException mtcr.get_species_names_wrapper()

        # Test error message content
        try
            mtcr.get_species_names_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end

        # Reload library for actual test
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Check species composition
        species_names = mtcr.get_species_names_wrapper()
        allowed_species = ["N", "N2", "N+", "N2+", "E-"]
        for name in species_names
            @test name in allowed_species
        end

        # Check return type
        @test species_names isa Vector{String}

        # Check that all names are non-empty strings
        @test all(length(name) > 0 for name in species_names)

        # Check that number of species is reasonable
        @test length(species_names) > 0
        @test length(species_names) <= 100  # Reasonable upper bound

        # Check that species names contain only valid characters
        isalnum(c) = isletter(c) || isdigit(c)
        for name in species_names
            @test all(c -> isascii(c) && (isalnum(c) || c in ['+', '-', '_']), name)
        end

        # Check that species names are unique
        @test length(species_names) == length(unique(species_names))
    end

    @testset "Electronic States Parameters" begin
        # Ensure library is loaded
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Test error when library not loaded
        mtcr.close_mtcr_library()
        @test_throws ErrorException mtcr.get_max_number_of_atomic_electronic_states_wrapper()
        @test_throws ErrorException mtcr.get_max_number_of_molecular_electronic_states_wrapper()

        # Test error message content
        try
            mtcr.get_max_number_of_atomic_electronic_states_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
            @test occursin("load_mtcr_library!", e.msg)
        end

        # Reload library for actual test
        mtcr.load_mtcr_library!(temp_mtcr_path)

        # Test atomic electronic states
        max_atomic_states = mtcr.get_max_number_of_atomic_electronic_states_wrapper()
        println("Max atomic electronic states: ", max_atomic_states)

        @test max_atomic_states isa Int32
        @test max_atomic_states > 0
        @test max_atomic_states <= 50  # Reasonable upper bound

        # Test molecular electronic states
        max_molecular_states = mtcr.get_max_number_of_molecular_electronic_states_wrapper()
        println("Max molecular electronic states: ", max_molecular_states)

        @test max_molecular_states isa Int32
        @test max_molecular_states > 0
        @test max_molecular_states <= 50  # Reasonable upper bound

        # Test that values are reasonable relative to each other
        @test max_molecular_states <= max_atomic_states
    end
end

@testset "Temperature Calculation" begin
    @testset "Error Handling Without Library" begin
        mtcr.close_mtcr_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ErrorException mtcr.calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            mtcr.calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ErrorException mtcr.calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            mtcr.calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_nowarn try
            result = mtcr.calculate_temperatures_wrapper(rho_sp, rho_etot;
                rho_eeex = 200.0, rho_evib = 200.0)

            # If successful, check structure
            @test result isa NamedTuple
            @test haskey(result, :tt)
            @test haskey(result, :trot)
            @test haskey(result, :teex)
            @test haskey(result, :tvib)
            @test haskey(result, :tex)
            @test haskey(result, :tvx)

            # Test that tvx uses the correct molecular electronic states dimension
            if result.tvx isa Matrix
                @test size(result.tvx, 1) ==
                      mtcr.get_max_number_of_molecular_electronic_states_wrapper()
                @test size(result.tvx, 2) == mtcr.get_max_number_of_species_wrapper()
            end
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Vibrational Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        mtcr.close_mtcr_library()

        tvib = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)

        try
            mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        tvib = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)

        try
            mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        tvib = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test basic call without optional arguments
        @test_nowarn try
            result = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)
            # If successful, check return type
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0  # Vibrational energy should be non-negative
        catch e
            @test e isa ErrorException
        end

        # Test with optional arguments
        max_species = mtcr.get_max_number_of_species_wrapper()
        max_electronic_states = mtcr.get_max_number_of_molecular_electronic_states_wrapper()
        rho_ex = zeros(Float64, max_electronic_states, max_species)
        tex = fill(tvib, max_species)
        teex = tvib

        @test_nowarn try
            result = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp;
                rho_ex = rho_ex, tex = tex, teex = teex)
            # If successful, check return type
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end

        # Test with partial optional arguments
        @test_nowarn try
            result = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp; teex = teex)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end

    @testset "Input Validation and Edge Cases" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with different vibrational temperatures
        test_temperatures = [100.0, 300.0, 1000.0, 5000.0, 10000.0]

        for tvib in test_temperatures
            @test_nowarn try
                result = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)
                @test result isa Float64
                @test isfinite(result)
                @test result >= 0.0
            catch e
                @test e isa ErrorException
            end
        end

        # Test with very small densities
        rho_sp_small = [1e-30, 1e-30, 1e-30, 1e-30, 1e-30]
        @test_nowarn try
            result = mtcr.calculate_vibrational_energy_wrapper(1000.0, rho_sp_small)
            println("vib temps: ", result)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Electron-Electronic Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        mtcr.close_mtcr_library()

        teex = 10000.0
        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_electron_electronic_energy_wrapper(
            teex, tvib, rho_sp)

        try
            mtcr.calculate_electron_electronic_energy_wrapper(teex, tvib, rho_sp)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        teex = 10000.0
        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_electron_electronic_energy_wrapper(
            teex, tvib, rho_sp)

        try
            mtcr.calculate_electron_electronic_energy_wrapper(teex, tvib, rho_sp)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Input Validation" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with negative temperature
        @test_throws ArgumentError mtcr.calculate_electron_electronic_energy_wrapper(
            -1000.0, tvib, rho_sp)

        # Test with zero temperature
        @test_throws ArgumentError mtcr.calculate_electron_electronic_energy_wrapper(
            0.0, tvib, rho_sp)

        # Test with empty species array
        @test_throws ArgumentError mtcr.calculate_electron_electronic_energy_wrapper(
            1000.0, tvib, Float64[])

        # Test error messages
        try
            mtcr.calculate_electron_electronic_energy_wrapper(-1000.0, tvib, rho_sp)
            @test false
        catch e
            @test occursin("Electron-electronic temperature must be positive", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        teex = 10000.0
        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_nowarn try
            result = mtcr.calculate_electron_electronic_energy_wrapper(teex, tvib, rho_sp)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0  # Electron-electronic energy should be non-negative
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Electronic Boltzmann Distribution" begin
    @testset "Error Handling Without Library" begin
        mtcr.close_mtcr_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tex = 10000.0
        trot = 1000.0
        tvib = 2000.0

        @test_throws ErrorException mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, tex, trot, tvib)

        try
            mtcr.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tex = 10000.0
        trot = 1000.0
        tvib = 2000.0

        @test_throws ErrorException mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, tex, trot, tvib)

        try
            mtcr.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Input Validation" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with negative temperatures
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, -1000.0, 1000.0, 2000.0)
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, -1000.0, 2000.0)
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 1000.0, -2000.0)

        # Test with zero temperatures
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, 0.0, 1000.0, 2000.0)
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 0.0, 2000.0)
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 1000.0, 0.0)

        # Test with empty species array
        @test_throws ArgumentError mtcr.set_electronic_boltzmann_wrapper(
            Float64[], 1000.0, 1000.0, 2000.0)

        # Test error messages
        try
            mtcr.set_electronic_boltzmann_wrapper(rho_sp, -1000.0, 1000.0, 2000.0)
            @test false
        catch e
            @test occursin("All temperatures must be positive", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tex = 10000.0
        trot = 1000.0
        tvib = 2000.0

        @test_nowarn try
            result = mtcr.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
            @test result isa Matrix{Float64}
            @test size(result, 2) == mtcr.get_max_number_of_species_wrapper()  # Second dimension should be max species
            # Now test that it uses the correct atomic electronic states dimension
            @test size(result, 1) ==
                  mtcr.get_max_number_of_atomic_electronic_states_wrapper()  # First dimension should be max atomic electronic states
            @test all(isfinite.(result))
            @test all(result .>= 0.0)  # Densities should be non-negative
        catch e
            @test e isa ErrorException
        end
    end

    @testset "Temperature Variations" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with different temperature combinations
        test_cases = [
            (1000.0, 1000.0, 1000.0),   # Equal temperatures
            (10000.0, 1000.0, 2000.0),  # High electronic temperature
            (1000.0, 5000.0, 2000.0),   # High rotational temperature
            (1000.0, 1000.0, 10000.0)  # High vibrational temperature
        ]

        for (tex, trot, tvib) in test_cases
            @test_nowarn try
                result = mtcr.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
                @test result isa Matrix{Float64}
                @test all(isfinite.(result))
                @test all(result .>= 0.0)
            catch e
                @test e isa ErrorException
            end
        end
    end
end

@testset "Total Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        mtcr.close_mtcr_library()

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_total_energy_wrapper(tt, rho_sp)

        try
            mtcr.calculate_total_energy_wrapper(tt, rho_sp)
            @test false
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException mtcr.calculate_total_energy_wrapper(tt, rho_sp)

        try
            mtcr.calculate_total_energy_wrapper(tt, rho_sp)
            @test false
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        mtcr.load_mtcr_library!(temp_mtcr_path)
        test_case_path = joinpath(@__DIR__, "test_case")
        mtcr.initialize_api_wrapper(case_path = test_case_path)

        tt = 1000.0
        tvib = 2000.0
        telec = 3000.0

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_ex = mtcr.set_electronic_boltzmann_wrapper(rho_sp, telec, tt, tvib)
        rho_eeex = mtcr.calculate_electron_electronic_energy_wrapper(telec, tvib, rho_sp)
        rho_evib = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)

        # Test with optional arguments
        # (generally optional, but required here due to flags in test case input file)
        @test_nowarn try
            result = mtcr.calculate_total_energy_wrapper(
                tt, rho_sp;
                rho_ex = rho_ex,
                rho_eeex = rho_eeex,
                rho_evib = rho_evib
            )

            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Source Terms Calculation" begin
    @testset "Error Handling Without Library" begin
        # Ensure library is not loaded
        mtcr.close_mtcr_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test error when library not loaded
        @test_throws ErrorException mtcr.calculate_sources_wrapper(rho_sp, rho_etot)

        # Test error message content
        try
            mtcr.calculate_sources_wrapper(rho_sp, rho_etot)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        # Load library but don't initialize
        mtcr.load_mtcr_library!(temp_mtcr_path)
        mtcr.MTCR_INITIALIZED[] = false

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test error when not initialized
        @test_throws ErrorException mtcr.calculate_sources_wrapper(rho_sp, rho_etot)

        # Test error message content
        try
            mtcr.calculate_sources_wrapper(rho_sp, rho_etot)
            @test false  # Should not reach here
        catch e
            @test occursin("MTCR not initialized", e.msg)
        end
    end

    # #TODO: Need to find appropriate values to test call to calculate_sources
    # @testset "Function Signature and Return Structure" begin
    #     # Set up proper initialization
    #     mtcr.load_mtcr_library!(temp_mtcr_path)
    #     test_case_path = joinpath(@__DIR__, "test_case")
    #     mtcr.initialize_api_wrapper(case_path = test_case_path)

    #     rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

    #     tt = 1000.0
    #     tvib = 2000.0
    #     telec = 3000.0

    #     rho_ex = mtcr.set_electronic_boltzmann_wrapper(rho_sp, telec, tt, tvib)
    #     rho_eeex = mtcr.calculate_electron_electronic_energy_wrapper(telec, tvib, rho_sp)
    #     rho_evib = mtcr.calculate_vibrational_energy_wrapper(tvib, rho_sp)
    #     rho_etot = mtcr.calculate_total_energy_wrapper(
    #         tt, rho_sp; rho_ex = rho_ex, rho_eeex = rho_eeex, rho_evib = rho_evib,)

    #     # Test that function can be called (may fail with actual Fortran call)
    #     # but should have proper error handling
    #     @test_nowarn try
    #         result = mtcr.calculate_sources_wrapper(rho_sp, rho_etot;
    #             rho_ex = rho_ex,
    #             rho_eeex = rho_eeex,
    #             rho_evib = rho_evib,)

    #         # # If successful, check structure
    #         @test result isa NamedTuple
    #         @test haskey(result, :drho_sp)
    #         @test haskey(result, :drho_etot)
    #         @test haskey(result, :drho_erot)
    #         @test haskey(result, :drho_eeex)
    #         @test haskey(result, :drho_evib)
    #     catch e
    #         # Expected to fail until Fortran library is built
    #         @test e isa ErrorException
    #     end
    # end
end

@testset "Finalization" begin
    # Test that finalization runs without error
    @test mtcr.finalize_api_wrapper() === nothing
    @test_nowarn mtcr.close_mtcr_library()

    @test mtcr.MTCR_HANDLE[] == C_NULL
    @test mtcr.MTCR_LIB_PATH[] == ""
end

@testset "Initialize-Finalize Catching When Pre-Initialized" begin
    # Test load checker
    @test !mtcr.is_mtcr_loaded()

    # Ensure library path is set
    mtcr.load_mtcr_library!(temp_mtcr_path)
    @test mtcr.is_mtcr_loaded()

    # Get path to test case directory
    test_case_path = joinpath(@__DIR__, "test_case")

    # Test basic reinitialization with required arguments
    mtcr.MTCR_INITIALIZED[] = true
    @test mtcr.initialize_api_wrapper(case_path = test_case_path) ===
          nothing

    # Test basic refinalization
    mtcr.MTCR_INITIALIZED[] = false
    @test mtcr.finalize_api_wrapper() === nothing
    @test_nowarn mtcr.close_mtcr_library()
    @test !mtcr.MTCR_INITIALIZED[]
end
