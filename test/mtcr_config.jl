
@testset "MTCR Configuration Tests" begin
    # Mock infrastructure for fortran_wrapper functions during testing
    mock_library_loaded = Ref{Bool}(false)
    mock_max_species = Ref{Int32}(20)
    mock_max_atomic_states = Ref{Int32}(10)
    mock_max_molecular_states = Ref{Int32}(5)
    mock_species_list = Ref{Vector{String}}(["N", "N2", "N+", "N2+", "E-"])

    # Helper functions to control mock behavior
    function set_mock_library_loaded(loaded::Bool)
        mock_library_loaded[] = loaded
    end

    function set_mock_mtcr_limits(max_species::Int, max_atomic::Int, max_molecular::Int)
        mock_max_species[] = Int32(max_species)
        mock_max_atomic_states[] = Int32(max_atomic)
        mock_max_molecular_states[] = Int32(max_molecular)
    end

    function set_mock_species_list(species::Vector{String})
        mock_species_list[] = species
    end

    @testset "TemperatureConfig" begin
        @testset "Valid Construction" begin
            # Test basic construction (keywords)
            temp_config = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            @test temp_config.Tt == 300.0
            @test temp_config.Tv == 300.0
            @test temp_config.Te == 10000.0
            @test temp_config.Tee == 300.0
        end

        @testset "Invalid Construction" begin
            # Test negative temperatures
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = -100.0, Tv = 300.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = -100.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = -400.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 400.0, Te = -10000.0)

            # Test zero temperatures
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 0.0, Tv = 300.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = 0.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 0.0, Te = 10000.0)
            @test_throws ErrorException mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 400.0, Te = 0.0)
        end
    end

    @testset "TimeIntegrationConfig" begin
        @testset "Valid Construction" begin
            # Test basic construction
            time_config = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = 2)
            @test time_config.dt == 1e-6
            @test time_config.dtm == 1e-4
            @test time_config.tlim == 1e-3
            @test time_config.nstep == 1000
            @test time_config.method == 2

            # Test with defaults
            time_config2 = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            @test time_config2.dt == 1e-6
            @test time_config2.dtm == 1e-4
            @test time_config2.tlim == 1e-3
            @test time_config2.nstep == 500000  # Default
            @test time_config2.method == 2      # Default
        end

        @testset "Invalid Construction" begin
            # Test negative time parameters
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = -1e-6, dtm = 1e-4, tlim = 1e-3)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = -1e-4, tlim = 1e-3)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = -1e-3)

            # Test zero time parameters
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 0.0, dtm = 1e-4, tlim = 1e-3)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 0.0, tlim = 1e-3)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 0.0)

            # Test invalid nstep
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 0)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = -100)

            # Test invalid method
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = 3)
            @test_throws ErrorException mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = -1)
        end
    end

    @testset "PhysicsConfig" begin
        @testset "Default Construction" begin
            physics = mtcr.PhysicsConfig()
            @test physics.bbh_model == 4
            @test physics.esc_model == 1
            @test physics.ar_et_model == 1
            @test physics.eex_noneq == 1
            @test physics.ev_relax_set == 1
            @test physics.et_relax_set == 1
        end

        @testset "Custom Construction" begin
            physics = mtcr.PhysicsConfig(
                bbh_model = 2,
                esc_model = 0,
                ar_et_model = 2,
                eex_noneq = 0,
                ev_relax_set = 2,
                et_relax_set = 2
            )
            @test physics.bbh_model == 2
            @test physics.esc_model == 0
            @test physics.ar_et_model == 2
            @test physics.eex_noneq == 0
            @test physics.ev_relax_set == 2
            @test physics.et_relax_set == 2
        end
    end

    @testset "ProcessConfig" begin
        @testset "Default Construction" begin
            processes = mtcr.ProcessConfig()
            @test processes.consider_elec_bbe == 1
            @test processes.consider_elec_bfe == 1
            @test processes.consider_elec_bbh == 1
            @test processes.consider_elec_bfh == 1
            @test processes.consider_rad == 0
            @test processes.consider_rdr == 0
            @test processes.consider_chem == 1
        end

        @testset "Custom Construction" begin
            processes = mtcr.ProcessConfig(
                consider_elec_bbe = 0,
                consider_elec_bfe = 0,
                consider_elec_bbh = 0,
                consider_elec_bfh = 0,
                consider_rad = 1,
                consider_rdr = 1,
                consider_chem = 0
            )
            @test processes.consider_elec_bbe == 0
            @test processes.consider_elec_bfe == 0
            @test processes.consider_elec_bbh == 0
            @test processes.consider_elec_bfh == 0
            @test processes.consider_rad == 1
            @test processes.consider_rdr == 1
            @test processes.consider_chem == 0
        end
    end

    @testset "MTCRConfig" begin
        @testset "Valid Construction" begin
            species = ["N", "N2", "N+", "N2+", "E-"]
            mole_fractions = [1e-20, 0.9998, 1e-20, 0.0001, 0.0001]
            total_number_density = 1e13
            temperatures = mtcr.TemperatureConfig(; Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 0.5e-5, dtm = 5.0, tlim = 1e3)

            config = mtcr.MTCRConfig(
                species = species,
                mole_fractions = mole_fractions,
                total_number_density = total_number_density,
                temperatures = temperatures,
                time_params = time_params
            )

            @test config.species == species
            @test config.mole_fractions == mole_fractions
            @test config.total_number_density == total_number_density
            @test config.temperatures == temperatures
            @test config.time_params == time_params
            @test config.unit_system == :CGS  # Default
            @test config.library_path == ""   # Default
            @test config.case_path == pwd()   # Default
            @test config.validate_species_against_mtcr == false  # Default
        end

        @testset "Custom Construction with All Parameters" begin
            species = ["Ar", "Ar+", "E-"]
            mole_fractions = [0.9998, 0.0001, 0.0001]
            total_number_density = 1e12
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 50000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-2)
            physics = mtcr.PhysicsConfig(bbh_model = 2)
            processes = mtcr.ProcessConfig(consider_rad = 1)
            test_case_path = joinpath(@__DIR__, "test_case")
            test_database_path = joinpath(test_case_path,
                "database/n2/elec_sts_expanded_electron_fits_ground")

            config = mtcr.MTCRConfig(
                species = species,
                mole_fractions = mole_fractions,
                total_number_density = total_number_density,
                temperatures = temperatures,
                time_params = time_params,
                physics = physics,
                processes = processes,
                database_path = test_database_path,
                library_path = "",
                case_path = test_case_path,
                unit_system = :SI,
                validate_species_against_mtcr = true,
                radiation_length = 2.0,
                print_source_terms = false,
                get_electron_density_by_charge_balance = false,
                min_sts_frac = 1e-25,
                is_isothermal_teex = false
            )

            @test config.species == species
            @test config.mole_fractions == mole_fractions
            @test config.total_number_density == total_number_density
            @test config.temperatures == temperatures
            @test config.time_params == time_params
            @test config.physics == physics
            @test config.processes == processes
            @test config.database_path == test_database_path
            @test config.library_path == ""
            @test config.case_path == test_case_path
            @test config.unit_system == :SI
            @test config.validate_species_against_mtcr == true
            @test config.radiation_length == 2.0
            @test config.print_source_terms == false
            @test config.get_electron_density_by_charge_balance == false
            @test config.min_sts_frac == 1e-25
            @test config.is_isothermal_teex == false
        end
    end

    @testset "validate_config" begin
        @testset "Valid Inputs" begin
            species = ["N", "N2", "E-"]
            mole_fractions = [0.1, 0.8, 0.1]
            total_number_density = 1e13
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            library_path = ""
            case_path = pwd()
            unit_system = :CGS

            @test mtcr.validate_config(species, mole_fractions, total_number_density,
                temperatures, time_params, library_path, case_path, unit_system) == true
        end

        @testset "Species and Mole Fractions Validation" begin
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            library_path = ""
            case_path = pwd()
            unit_system = :CGS

            # Test mismatched array lengths
            @test_throws ArgumentError mtcr.validate_config(
                ["N", "N2"], [0.5], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test empty species
            @test_throws ArgumentError mtcr.validate_config(
                String[], Float64[], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test mole fractions don't sum to 1
            @test_throws ArgumentError mtcr.validate_config(
                ["N", "N2"], [0.3, 0.3], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test negative mole fractions
            @test_throws ArgumentError mtcr.validate_config(
                ["N", "N2"], [-0.1, 1.1], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test duplicate species
            @test_throws ArgumentError mtcr.validate_config(
                ["N", "N"], [0.5, 0.5], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test empty species name
            @test_throws ArgumentError mtcr.validate_config(
                ["N", ""], [0.5, 0.5], 1e13, temperatures, time_params,
                library_path, case_path, unit_system)
        end

        @testset "Number Density Validation" begin
            species = ["N", "N2"]
            mole_fractions = [0.5, 0.5]
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            library_path = ""
            case_path = pwd()
            unit_system = :CGS

            # Test negative number density
            @test_throws ArgumentError mtcr.validate_config(
                species, mole_fractions, -1e13, temperatures, time_params,
                library_path, case_path, unit_system)

            # Test zero number density
            @test_throws ArgumentError mtcr.validate_config(
                species, mole_fractions, 0.0, temperatures, time_params,
                library_path, case_path, unit_system)
        end

        @testset "Path Validation" begin
            species = ["N", "N2"]
            mole_fractions = [0.5, 0.5]
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            unit_system = :CGS

            # Test non-existent library path
            @test_throws ArgumentError mtcr.validate_config(
                species, mole_fractions, 1e13, temperatures, time_params,
                "/nonexistent/library.so", pwd(), unit_system)

            # Test non-existent case path
            @test_throws ArgumentError mtcr.validate_config(
                species, mole_fractions, 1e13, temperatures, time_params,
                "", "/nonexistent/path", unit_system)
        end

        @testset "Unit System Validation" begin
            species = ["N", "N2"]
            mole_fractions = [0.5, 0.5]
            temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            library_path = ""
            case_path = pwd()

            # Test invalid unit system
            @test_throws ArgumentError mtcr.validate_config(
                species, mole_fractions, 1e13, temperatures, time_params,
                library_path, case_path, :INVALID)
        end
    end

    @testset "nitrogen_10ev_config" begin
        @testset "Default Configuration" begin
            config = mtcr.nitrogen_10ev_config()

            @test config.species == ["N", "N2", "N+", "N2+", "E-"]
            @test config.mole_fractions == [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
            @test config.total_number_density == 1.0e13
            @test config.temperatures.Tt == 750.0
            @test config.temperatures.Tv == 750.0
            @test config.temperatures.Te == 115000.0
            @test config.temperatures.Tee == 750.0
            # Stored in seconds (prob_setup writes microseconds)
            @test config.time_params.dt ≈ 5e-12
            @test config.time_params.dtm ≈ 5e-6
            @test config.time_params.tlim ≈ 1e-3
            @test config.time_params.nstep == 500000
            @test config.time_params.method == 2
        end
    end

    @testset "generate_input_files" begin
        @testset "Directory Structure Creation" begin
            # Create a temporary directory for testing
            temp_dir = mktempdir()
            try
                config = mtcr.nitrogen_10ev_config()

                # Test successful file generation
                @test mtcr.generate_input_files(config, temp_dir) == true

                # Check that required directories were created
                @test isdir(joinpath(temp_dir, "input"))
                @test isdir(joinpath(temp_dir, "output"))
                @test isdir(joinpath(temp_dir, "output", "sources"))
                @test isdir(joinpath(temp_dir, "output", "states"))

                # Check that input files were created
                @test isfile(joinpath(temp_dir, "input", "prob_setup.inp"))
                @test isfile(joinpath(temp_dir, "input", "sources_setup.inp"))
                @test isfile(joinpath(temp_dir, "input", "tau_scaling.inp"))

            finally
                rm(temp_dir; recursive = true)
            end
        end

        @testset "File Content Validation" begin
            temp_dir = mktempdir()
            try
                config = mtcr.nitrogen_10ev_config()
                mtcr.generate_input_files(config, temp_dir)

                # Read and check prob_setup.inp content
                prob_setup_content = read(
                    joinpath(temp_dir, "input", "prob_setup.inp"), String)
                @test occursin("NSP=5", prob_setup_content)  # 5 species

                @test occursin("X1=1.0e-20", prob_setup_content)
                @test occursin("X2=0.9998", prob_setup_content)
                @test occursin("X3=1.0e-20", prob_setup_content)
                @test occursin("X4=0.0001", prob_setup_content)
                @test occursin("X5=0.0001", prob_setup_content)

                @test occursin("TOTAL_NUMBER_DENSITY=1.0e13", prob_setup_content)

                @test occursin("TT=750.0", prob_setup_content)
                @test occursin("TV=750.0", prob_setup_content)
                @test occursin("TEE=750.0", prob_setup_content)
                @test occursin("TE=115000.0", prob_setup_content)

                @test occursin("RAD_LEN=1.0", prob_setup_content)

                @test occursin("BBHMODEL=4", prob_setup_content)
                @test occursin("ESC_MODEL=1", prob_setup_content)
                @test occursin("AR_ET_MODEL=1", prob_setup_content)
                @test occursin("EEX_NONEQ=1", prob_setup_content)
                @test occursin("EV_RELAX_SET=1", prob_setup_content)
                @test occursin("ET_RELAX_SET=1", prob_setup_content)

                @test occursin("CONSIDER_ELEC_BBE=1", prob_setup_content)
                @test occursin("CONSIDER_ELEC_BFE=1", prob_setup_content)
                @test occursin("CONSIDER_ELEC_BBH=1", prob_setup_content)
                @test occursin("CONSIDER_ELEC_BFH=1", prob_setup_content)
                @test occursin("CONSIDER_RAD=0", prob_setup_content)
                @test occursin("CONSIDER_RDR=0", prob_setup_content)
                @test occursin("CONSIDER_CHEM=1", prob_setup_content)

                @test occursin("TIME_METHOD=2", prob_setup_content)
                @test occursin("IS_ISOTHERMAL_TEEX=0", prob_setup_content)
                @test occursin("ND=0", prob_setup_content)
                @test occursin("DT=5.0e-6", prob_setup_content)
                @test occursin("DTM=5.0", prob_setup_content)
                @test occursin("TLIM=1000.0", prob_setup_content)
                @test occursin("NSTEP=500000", prob_setup_content)

                # Read and check sources_setup.inp content
                sources_content = read(
                    joinpath(temp_dir, "input", "sources_setup.inp"), String)
                @test occursin("BEGIN SPECIES SOURCES", sources_content)
                @test occursin("END SPECIES SOURCES", sources_content)
                @test occursin("BEGIN EXCITED STATE SOURCES", sources_content)
                @test occursin("END EXCITED STATE SOURCES", sources_content)
                for species in config.species
                    @test occursin(species, sources_content)
                end

            finally
                rm(temp_dir; recursive = true)
            end
        end
    end

    @testset "get_molecular_weights" begin
        @testset "Known Species" begin
            # Test nitrogen species
            weights = mtcr.get_molecular_weights(["N", "N2", "N+", "N2+", "E-"])
            @test weights[1] ≈ 14.007  # N
            @test weights[2] ≈ 28.014  # N2
            @test weights[3] ≈ 14.007  # N+
            @test weights[4] ≈ 28.014  # N2+
            @test weights[5] ≈ 5.485799e-4  # E-

            # Test noble gas species
            weights_ar = mtcr.get_molecular_weights(["Ar", "Ar+"])
            @test weights_ar[1] ≈ 39.948  # Ar
            @test weights_ar[2] ≈ 39.948  # Ar+

            weights_xe = mtcr.get_molecular_weights(["Xe", "Xe+"])
            @test weights_xe[1] ≈ 131.293  # Xe
            @test weights_xe[2] ≈ 131.293  # Xe+
        end

        @testset "Unknown Species" begin
            @test_throws ErrorException mtcr.get_molecular_weights(["UNKNOWN"])
            @test_throws ErrorException mtcr.get_molecular_weights(["N", "UNKNOWN", "E-"])
        end
    end

    @testset "MTCR Fortran Interface Validation Tests" begin
        @testset "validate_config_against_mtcr" begin
            @testset "Basic Functionality" begin
                config = mtcr.MTCRConfig(
                    species = ["N", "N2", "E-"],
                    mole_fractions = [0.1, 0.8, 0.1],
                    total_number_density = 1e13,
                    temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                    time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
                )

                # Should return true (may show warnings if library not loaded)
                @test mtcr.validate_config_against_mtcr(config) == true
            end

            @testset "With Species Validation Enabled" begin
                config = mtcr.MTCRConfig(
                    species = ["N", "N2", "E-"],
                    mole_fractions = [0.1, 0.8, 0.1],
                    total_number_density = 1e13,
                    temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                    time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                    validate_species_against_mtcr = true
                )

                # Should return true (may show warnings if library not loaded or species not found)
                @test mtcr.validate_config_against_mtcr(config) == true
            end
        end

        @testset "validate_species_against_mtcr_database (Enhanced)" begin
            @testset "Validation Disabled" begin
                config = mtcr.MTCRConfig(
                    species = ["N", "N2", "E-"],
                    mole_fractions = [0.1, 0.8, 0.1],
                    total_number_density = 1e13,
                    temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                    time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                    validate_species_against_mtcr = false
                )

                @test mtcr.validate_species_against_mtcr_database(config) == true
            end

            @testset "Validation Enabled" begin
                config = mtcr.MTCRConfig(
                    species = ["N", "N2", "E-"],
                    mole_fractions = [0.1, 0.8, 0.1],
                    total_number_density = 1e13,
                    temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                    time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                    validate_species_against_mtcr = true
                )

                # Should return true (may show warnings if library not loaded or species not found)
                @test mtcr.validate_species_against_mtcr_database(config) == true
            end
        end
    end

    @testset "convert_config_units" begin
        @testset "No Conversion Needed" begin
            config = mtcr.nitrogen_10ev_config()  # Default is CGS
            converted = mtcr.convert_config_units(config, :CGS)
            @test converted === config  # Should return same object
        end

        @testset "SI to CGS Conversion" begin
            # Create SI config
            config_si = mtcr.MTCRConfig(
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e19,  # SI: 1/m³
                temperatures = mtcr.TemperatureConfig(; Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                time_params = mtcr.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                unit_system = :SI
            )

            converted = mtcr.convert_config_units(config_si, :CGS)
            @test converted.unit_system == :CGS
            @test converted.total_number_density ≈ 1e13  # Converted to 1/cm³
            @test converted.species == config_si.species
            @test converted.mole_fractions == config_si.mole_fractions
        end

        @testset "CGS to SI Conversion" begin
            config_cgs = mtcr.nitrogen_10ev_config()  # CGS by default
            converted = mtcr.convert_config_units(config_cgs, :SI)

            @test converted.unit_system == :SI
            @test converted.total_number_density ≈ 1e19  # Converted to 1/m³
            @test converted.species == config_cgs.species
            @test converted.mole_fractions == config_cgs.mole_fractions
        end

        @testset "Invalid Conversion" begin
            config = mtcr.nitrogen_10ev_config()
            @test_throws ErrorException mtcr.convert_config_units(config, :INVALID)
        end
    end

    @testset "MTCRResults" begin
        @testset "Structure Creation" begin
            time = [0.0, 1.0, 2.0]
            species_densities = [1e-3 1e-3 1e-3; 1e-6 1e-6 1e-6]
            temperatures = (tt = [300.0, 310.0, 320.0], te = [10000.0, 11000.0, 12000.0])
            total_energy = [1e4, 1.1e4, 1.2e4]

            results = mtcr.MTCRResults(
                time, species_densities, temperatures, total_energy,
                nothing, true, "Success"
            )

            @test results.time == time
            @test results.species_densities == species_densities
            @test results.temperatures == temperatures
            @test results.total_energy == total_energy
            @test results.source_terms === nothing
            @test results.success == true
            @test results.message == "Success"
        end
    end
end
