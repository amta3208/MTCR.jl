"""
# Data Conversion Module

This module handles conversion between Julia and Fortran data structures,
including unit conversions and array format transformations.

## Key Functions

- `julia_to_fortran_arrays`: Convert Julia arrays to Fortran-compatible format
- `fortran_to_julia_arrays`: Convert Fortran arrays back to Julia format
- `convert_units_si_to_cgs`: Convert from SI to CGS units (MTCR uses CGS)
- `convert_units_cgs_to_si`: Convert from CGS to SI units
- `validate_species_data`: Validate species data consistency

## Unit Conversions

MTCR uses CGS units while HallThruster.jl typically uses SI units.
This module handles all necessary conversions:

- Density: kg/m³ ↔ g/cm³
- Energy: J ↔ erg
- Temperature: K (same in both systems)
- Pressure: Pa ↔ dyne/cm²
- Time: s (same in both systems)
"""

# Unit conversion constants
const KG_M3_TO_G_CM3 = 1e-3  # kg/m³ to g/cm³
const G_CM3_TO_KG_M3 = 1e3   # g/cm³ to kg/m³
const J_TO_ERG = 1e7         # J to erg
const ERG_TO_J = 1e-7        # erg to J
const PA_TO_DYNE_CM2 = 10.0  # Pa to dyne/cm²
const DYNE_CM2_TO_PA = 0.1   # dyne/cm² to Pa
const J_M3_TO_ERG_CM3 = 10.0  # J/m³ to erg/cm³
const ERG_CM3_TO_J_M3 = 0.1   # erg/cm³ to J/m³
const AVOGADRO = 6.02214076e23  # Avogadro's number (1/mol)

"""
$(SIGNATURES)

Convert species densities from SI (kg/m³) to CGS (g/cm³) units.
"""
function convert_density_si_to_cgs(rho_si::Vector{Float64})
    return rho_si .* KG_M3_TO_G_CM3
end

"""
$(SIGNATURES)

Convert species densities from CGS (g/cm³) to SI (kg/m³) units.
"""
function convert_density_cgs_to_si(rho_cgs::Vector{Float64})
    return rho_cgs .* G_CM3_TO_KG_M3
end

"""
$(SIGNATURES)

Convert energy density from SI (J/m³) to CGS (erg/cm³) units.
"""
function convert_energy_density_si_to_cgs(energy_si::Float64)
    return energy_si * J_M3_TO_ERG_CM3
end

"""
$(SIGNATURES)

Convert energy density from CGS (erg/cm³) to SI (J/m³) units.
"""
function convert_energy_density_cgs_to_si(energy_cgs::Float64)
    return energy_cgs * ERG_CM3_TO_J_M3
end

"""
$(SIGNATURES)

Convert pressure from SI (Pa) to CGS (dyne/cm²) units.
"""
function convert_pressure_si_to_cgs(pressure_si::Float64)
    return pressure_si * PA_TO_DYNE_CM2
end

"""
$(SIGNATURES)

Convert pressure from CGS (dyne/cm²) to SI (Pa) units.
"""
function convert_pressure_cgs_to_si(pressure_cgs::Float64)
    return pressure_cgs * DYNE_CM2_TO_PA
end

"""
$(SIGNATURES)

Convert number density from SI (1/m³) to CGS (1/cm³) units.
"""
function convert_number_density_si_to_cgs(n_si::Float64)
    return n_si * 1e-6  # 1/m³ to 1/cm³
end

"""
$(SIGNATURES)

Convert number density from CGS (1/cm³) to SI (1/m³) units.
"""
function convert_number_density_cgs_to_si(n_cgs::Float64)
    return n_cgs * 1e6  # 1/cm³ to 1/m³
end

"""
$(SIGNATURES)

Convert a complete state vector from SI to CGS units for MTCR input.

# Arguments
- `rho_sp_si::Vector{Float64}`: Species densities in SI units (kg/m³)
- `rho_etot_si::Float64`: Total energy density in SI units (J/m³)
- `number_density_si::Float64`: Total number density in SI units (1/m³)
- Additional optional energy components in SI units

# Returns
- Named tuple with all quantities converted to CGS units
"""
function convert_state_si_to_cgs(rho_sp_si::Vector{Float64},
        rho_etot_si::Float64,
        number_density_si::Float64;
        rho_erot_si::Union{Float64, Nothing} = nothing,
        rho_eeex_si::Union{Float64, Nothing} = nothing,
        rho_evib_si::Union{Float64, Nothing} = nothing)
    rho_sp_cgs = convert_density_si_to_cgs(rho_sp_si)
    rho_etot_cgs = convert_energy_density_si_to_cgs(rho_etot_si)
    number_density_cgs = convert_number_density_si_to_cgs(number_density_si)

    rho_erot_cgs = rho_erot_si !== nothing ? convert_energy_density_si_to_cgs(rho_erot_si) :
                   nothing
    rho_eeex_cgs = rho_eeex_si !== nothing ? convert_energy_density_si_to_cgs(rho_eeex_si) :
                   nothing
    rho_evib_cgs = rho_evib_si !== nothing ? convert_energy_density_si_to_cgs(rho_evib_si) :
                   nothing

    return (rho_sp = rho_sp_cgs,
        rho_etot = rho_etot_cgs,
        number_density = number_density_cgs,
        rho_erot = rho_erot_cgs,
        rho_eeex = rho_eeex_cgs,
        rho_evib = rho_evib_cgs)
end

"""
$(SIGNATURES)

Convert a complete state vector from CGS to SI units from MTCR output.

# Arguments
- `rho_sp_cgs::Vector{Float64}`: Species densities in CGS units (g/cm³)
- `rho_etot_cgs::Float64`: Total energy density in CGS units (erg/cm³)
- Additional optional energy components in CGS units

# Returns
- Named tuple with all quantities converted to SI units
"""
function convert_state_cgs_to_si(rho_sp_cgs::Vector{Float64},
        rho_etot_cgs::Float64;
        rho_erot_cgs::Union{Float64, Nothing} = nothing,
        rho_eeex_cgs::Union{Float64, Nothing} = nothing,
        rho_evib_cgs::Union{Float64, Nothing} = nothing)
    rho_sp_si = convert_density_cgs_to_si(rho_sp_cgs)
    rho_etot_si = convert_energy_density_cgs_to_si(rho_etot_cgs)

    rho_erot_si = rho_erot_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_erot_cgs) : nothing
    rho_eeex_si = rho_eeex_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_eeex_cgs) : nothing
    rho_evib_si = rho_evib_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_evib_cgs) : nothing

    return (rho_sp = rho_sp_si,
        rho_etot = rho_etot_si,
        rho_erot = rho_erot_si,
        rho_eeex = rho_eeex_si,
        rho_evib = rho_evib_si)
end

"""
$(SIGNATURES)

Convert source terms from CGS to SI units.

# Arguments
- `drho_sp_cgs::Vector{Float64}`: Species source terms in CGS units (g/cm³/s)
- `drho_etot_cgs::Float64`: Energy source term in CGS units (erg/cm³/s)
- Additional optional source terms in CGS units

# Returns
- Named tuple with all source terms converted to SI units
"""
function convert_sources_cgs_to_si(drho_sp_cgs::Vector{Float64},
        drho_etot_cgs::Float64;
        drho_erot_cgs::Union{Float64, Nothing} = nothing,
        drho_eeex_cgs::Union{Float64, Nothing} = nothing,
        drho_evib_cgs::Union{Float64, Nothing} = nothing)

    # Source terms have units of [quantity]/time, so same conversion as state
    drho_sp_si = convert_density_cgs_to_si(drho_sp_cgs)
    drho_etot_si = convert_energy_density_cgs_to_si(drho_etot_cgs)

    drho_erot_si = drho_erot_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_erot_cgs) : nothing
    drho_eeex_si = drho_eeex_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_eeex_cgs) : nothing
    drho_evib_si = drho_evib_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_evib_cgs) : nothing

    return (drho_sp = drho_sp_si,
        drho_etot = drho_etot_si,
        drho_erot = drho_erot_si,
        drho_eeex = drho_eeex_si,
        drho_evib = drho_evib_si)
end

"""
$(SIGNATURES)

Prepare Julia arrays for passing to Fortran.

Fortran expects column-major arrays, which Julia uses natively.
This function ensures proper memory layout and type consistency.

# Arguments
- `arrays...`: Variable number of Julia arrays to prepare

# Returns
- Tuple of arrays ready for Fortran ccall
"""
function prepare_arrays_for_fortran(arrays...)
    prepared = []

    for arr in arrays
        if arr === nothing
            push!(prepared, nothing)
        elseif isa(arr, AbstractArray)
            # Ensure contiguous memory layout and Float64 type
            prepared_arr = Array{Float64}(arr)
            push!(prepared, prepared_arr)
        else
            # Scalar values
            push!(prepared, Float64(arr))
        end
    end

    return tuple(prepared...)
end

"""
$(SIGNATURES)

Validate species data consistency between Julia and MTCR.

# Arguments
- `species_names::Vector{String}`: Species names from configuration
- `mtcr_species::Vector{String}`: Species names from MTCR
- `densities::Vector{Float64}`: Species densities

# Returns
- `true` if validation passes, throws error otherwise
"""
function validate_species_data(species_names::Vector{String},
        mtcr_species::Vector{String},
        densities::Vector{Float64})

    # Check array lengths match
    if length(species_names) != length(densities)
        error("Species names and densities arrays have different lengths: $(length(species_names)) vs $(length(densities))")
    end

    # Check for negative densities
    if any(densities .< 0)
        negative_indices = findall(densities .< 0)
        error("Negative densities found for species: $(species_names[negative_indices])")
    end

    # Check species names are valid
    for (i, name) in enumerate(species_names)
        if !(name in mtcr_species)
            error("Species '$(name)' not found in MTCR database. Available species: $(mtcr_species)")
        end
    end

    # Check for very small densities that might cause numerical issues
    min_density = 1e-30
    small_density_indices = findall(densities .< min_density)
    if !isempty(small_density_indices)
        @warn "Very small densities detected for species: $(species_names[small_density_indices]). This may cause numerical issues."
    end

    return true
end

"""
$(SIGNATURES)

Create species mapping between HallThruster.jl and MTCR conventions.

# Arguments
- `ht_species::Vector{String}`: Species names in HallThruster.jl format
- `mtcr_species::Vector{String}`: Species names in MTCR format

# Returns
- Dictionary mapping HallThruster.jl species names to MTCR species names
"""
function create_species_mapping(ht_species::Vector{String}, mtcr_species::Vector{String})
    mapping = Dict{String, String}()

    # Common mappings (this may need to be expanded based on actual usage)
    common_mappings = Dict(
        "N" => "N",
        "N2" => "N2",
        "N+" => "N+",
        "N2+" => "N2+",
        "e-" => "E-",
        "E-" => "E-",
        "Ar" => "Ar",
        "Ar+" => "Ar+",
        "Xe" => "Xe",
        "Xe+" => "Xe+",
        "Kr" => "Kr",
        "Kr+" => "Kr+"
    )

    for ht_name in ht_species
        if haskey(common_mappings, ht_name)
            mtcr_name = common_mappings[ht_name]
            if mtcr_name in mtcr_species
                mapping[ht_name] = mtcr_name
            else
                error("MTCR species '$(mtcr_name)' not found in database for HallThruster species '$(ht_name)'")
            end
        else
            # Try direct mapping
            if ht_name in mtcr_species
                mapping[ht_name] = ht_name
            else
                error("No mapping found for HallThruster species '$(ht_name)' to MTCR database")
            end
        end
    end

    return mapping
end

"""
$(SIGNATURES)

Convert mole fractions to mass densities.

# Arguments
- `mole_fractions::Vector{Float64}`: Species mole fractions
- `molecular_weights::Vector{Float64}`: Species molecular weights (g/mol)
- `total_number_density::Float64`: Total number density (1/cm³)

# Returns
- Vector of mass densities (g/cm³)
"""
function mole_fractions_to_mass_densities(mole_fractions::Vector{Float64},
        molecular_weights::Vector{Float64},
        total_number_density::Float64)
    if length(mole_fractions) != length(molecular_weights)
        error("Mole fractions and molecular weights arrays must have same length")
    end

    if abs(sum(mole_fractions) - 1.0) > 1e-10
        error("Mole fractions must sum to 1.0, got $(sum(mole_fractions))")
    end

    # Convert to number densities
    number_densities = mole_fractions .* total_number_density

    # Convert to mass densities (g/cm³)
    # number_density (1/cm³) * molecular_weight (g/mol) / Avogadro (1/mol) = mass_density (g/cm³)
    mass_densities = number_densities .* molecular_weights ./ AVOGADRO

    return mass_densities
end

"""
$(SIGNATURES)

Convert mass densities to mole fractions.

# Arguments
- `mass_densities::Vector{Float64}`: Species mass densities (g/cm³)
- `molecular_weights::Vector{Float64}`: Species molecular weights (g/mol)

# Returns
- Vector of mole fractions
"""
function mass_densities_to_mole_fractions(mass_densities::Vector{Float64},
        molecular_weights::Vector{Float64})
    if length(mass_densities) != length(molecular_weights)
        error("Mass densities and molecular weights arrays must have same length")
    end

    # Convert to number densities
    number_densities = mass_densities .* AVOGADRO ./ molecular_weights

    # Convert to mole fractions
    total_number_density = sum(number_densities)

    if total_number_density ≈ 0.0
        error("Total number density is zero - cannot compute mole fractions")
    end

    mole_fractions = number_densities ./ total_number_density

    return mole_fractions
end
