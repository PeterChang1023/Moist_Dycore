module Restart_Manager_Module

using JLD2
using ..Dyn_Data_Module
using ..Experiment_Configuration
using ..Vert_Coordinate_Module
using ..Atmo_Data_Module

export Restart_Manager, Write_Restart_File, Load_Restart_File!, Cleanup_Old_Restarts, Restart_Contract

const RESTART_FORMAT_VERSION = 4

const REQUIRED_RESTART_FIELDS = (
    :spe_vor_p, :spe_vor_c, :spe_vor_n, :spe_div_p, :spe_div_c, :spe_div_n,
    :spe_lnps_p, :spe_lnps_c, :spe_lnps_n, :spe_t_p, :spe_t_c, :spe_t_n,
    :grid_u_p, :grid_u_c, :grid_u_n, :grid_v_p, :grid_v_c, :grid_v_n,
    :grid_ps_p, :grid_ps_c, :grid_ps_n, :grid_t_p, :grid_t_c, :grid_t_n,
    :grid_q_p, :grid_q_c, :grid_q_n, :grid_geopots,
)

"""Configuration needed to interpret a saved leapfrog state.

Forcing choices (BM, LRF, SST, etc.) and output settings may intentionally change
on restart; equation set, time history, grid and core dynamics must match.
"""
function Restart_Contract(config::Model_Config, vert_coord, atmo_data::Atmo_Data)
    contract = Dict{String,Any}("contract_version" => 1, "integrator" => "Filtered_Leapfrog")
    for key in (:model_type, :num_fourier, :nθ, :nd, :Δt, :day_to_sec,
                :damping_order, :damping_coef, :robert_coef, :implicit_coef,
                :moisture_processes)
        contract[string(key)] = getfield(config, key)
    end
    for key in (:radius, :omega, :grav, :rdgas, :kappa, :rvgas, :cp_air, :Lv,
                :gamma, :use_virtual_temperature, :do_mass_correction,
                :do_energy_correction, :do_water_correction)
        contract[string(key)] = getfield(atmo_data, key)
    end
    if vert_coord !== nothing
        for key in (:vert_coord_option, :vert_difference_option, :vert_advect_scheme,
                    :p_ref, :zero_top, :ak, :bk)
            contract[string(key)] = deepcopy(getfield(vert_coord, key))
        end
    end
    return contract
end

"""Construct a contract from the original configuration for attested legacy files."""
function Restart_Contract(config::Model_Config)
    nλ = 2config.nθ
    vert = config.model_type == :PrimitiveEquation ? Vert_Coordinate(nλ,config.nθ,config.nd,
        config.vert_coord_option,config.vert_difference_option,config.vert_ref_level_option) : nothing
    p = config.physics_params
    atmo = Atmo_Data(config.name,nλ,config.nθ,config.nd,
        get(p,"do_mass_correction",true),get(p,"do_energy_correction",true),
        get(p,"do_water_correction",true),get(p,"use_virtual_temperature",true),
        zeros(config.nθ); radius=config.radius,omega=config.omega,grav=config.grav,
        [Symbol(k)=>v for (k,v) in p]...)
    Restart_Contract(config,vert,atmo)
end


struct Restart_Manager
    output_dir::String
    saving_frequency::Int64
    contract::Union{Nothing,Dict{String,Any}}

    function Restart_Manager(output_dir::String, frequency::Int64; contract=nothing)
        if frequency > 0 && !isdir(output_dir)
            mkpath(output_dir)
        end
        return new(output_dir, frequency, deepcopy(contract))
    end
end



"""
    Write_Restart_File(manager, dyn_data, current_time)
"""
function Write_Restart_File(
    manager::Restart_Manager,
    dyn_data::Dyn_Data,
    current_time::Int64,
)
    # Ensure directory exists right before writing (lazy creation)
    if !isdir(manager.output_dir)
        mkpath(manager.output_dir)
    end

    filename = joinpath(manager.output_dir, "restart_t$(current_time).jld2")
    temp_filename = filename * ".tmp"

    jldopen(temp_filename, "w") do file
        # The unconfigured low-level writer retains v3 for raw state inspection.
        # Driver checkpoints always carry a v4 compatibility contract.
        file["restart_format_version"] = manager.contract === nothing ? 3 : RESTART_FORMAT_VERSION
        manager.contract === nothing || (file["restart_contract"] = manager.contract)
        file["saved_time"] = current_time
        file["dimensions"] =
            (dyn_data.num_fourier, dyn_data.num_spherical, dyn_data.nλ, dyn_data.nθ, dyn_data.nd)
        for name in fieldnames(Dyn_Data)
            value = getfield(dyn_data, name)
            value isa Array && (file["state/$(name)"] = value)
        end
    end

    mv(temp_filename, filename; force = true)
    @info "Checkpoint saved: $filename"
end



"""
    Load_Restart_File!(dyn_data, filename; expected_contract, legacy_contract, allow_partial)

Validate the whole state before copying. `allow_partial=true` is exclusively
for low-level inspection of incomplete historical files, never a forecast
restart. High-level initialization always supplies `expected_contract`.
"""
function Load_Restart_File!(dyn_data::Dyn_Data, filename::String;
                            expected_contract=nothing, legacy_contract=nothing,
                            allow_partial::Bool=false)
    allow_partial && expected_contract !== nothing &&
        error("Partial state loading cannot be used for a configured forecast restart")
    if !isfile(filename)
        error("Restart file not found: $filename")
    end

    @info "Loading warm start from: $filename"

    return jldopen(filename, "r") do file
        saved_time = file["saved_time"]
        saved_time isa Integer && !(saved_time isa Bool) &&
            0 <= saved_time <= typemax(Int64) || error("Invalid checkpoint saved_time")
        if expected_contract !== nothing
            contract = haskey(file,"restart_contract") ? file["restart_contract"] : legacy_contract
            contract === nothing && error("Checkpoint has no compatibility contract; " *
                "supply restart_legacy_contract = Restart_Contract(original_config) only after " *
                "verifying the original run configuration")
            contract isa AbstractDict || error("Invalid checkpoint compatibility contract")
            for (key,value) in expected_contract
                haskey(contract,key) && isequal(contract[key],value) ||
                    error("Restart configuration mismatch for $key")
            end
            dt = expected_contract["Δt"]
            saved_time >= dt && saved_time % dt == 0 ||
                error("Restart timestamp is incompatible with mature leapfrog history")
        end
        arrays = Dict{Symbol,Any}()
        if haskey(file, "restart_format_version")
            version = Int(file["restart_format_version"])
            version in (2, 3, RESTART_FORMAT_VERSION) ||
                error("Unsupported restart format version $version")
            version == 4 && !haskey(file,"restart_contract") &&
                error("Version-4 checkpoint is missing its compatibility contract")
            expected =
                (dyn_data.num_fourier, dyn_data.num_spherical, dyn_data.nλ, dyn_data.nθ, dyn_data.nd)
            file["dimensions"] == expected ||
                error("Restart resolution does not match the configured model")
            for name in fieldnames(Dyn_Data)
                path = "state/$(name)"
                haskey(file, path) || continue
                getfield(dyn_data,name) isa Array || continue
                arrays[name] = file[path]
            end
        else
            # Legacy files stored the complete Dyn_Data object. Copy common
            # array fields so removed scratch fields are harmless.
            legacy = file["dyn_data_state"]
            for name in fieldnames(Dyn_Data)
                hasproperty(legacy, name) || continue
                dest = getfield(dyn_data, name)
                dest isa Array || continue
                arrays[name] = getproperty(legacy, name)
            end
        end
        if !allow_partial
            missing = filter(name -> !haskey(arrays,name), REQUIRED_RESTART_FIELDS)
            isempty(missing) || error("Checkpoint missing required state/history fields: $missing")
        end
        # Validation precedes every copy, so an invalid late field leaves the
        # destination completely untouched and preserves all existing buffers.
        for (name,src) in arrays
            dest = getfield(dyn_data,name)
            src isa Array || error("Restart field $name is not an array")
            size(src) == size(dest) || error("Restart field $name has size $(size(src)); expected $(size(dest))")
            eltype(src) == eltype(dest) || error("Restart field $name has incompatible element type")
            all(isfinite,src) || error("Restart field $name contains nonfinite values")
        end
        if expected_contract !== nothing && expected_contract["model_type"] == :PrimitiveEquation
            for level in (:p,:c,:n), stem in (:grid_ps,:grid_t)
                name = Symbol(stem,"_",level)
                all(>(0),arrays[name]) || error("Restart field $name must be positive")
            end
            for level in (:p,:c,:n)
                name = Symbol("grid_q_",level)
                all(q -> -1e-12 <= q <= 1,arrays[name]) ||
                    error("Restart field $name has invalid specific humidity")
            end
        end
        for (name,src) in arrays
            copyto!(getfield(dyn_data,name),src)
        end
        Int64(saved_time)
    end
end


"""
    Cleanup_Old_Restarts(manager, keep_last_n, protect_file)
    
Deletes older restart files, BUT ensures the file used to start the run 
(protect_file) is never deleted.
"""
function Cleanup_Old_Restarts(
    manager::Restart_Manager,
    keep_last_n::Int = 3,
    protect_file::String = "",
)
    files = readdir(manager.output_dir, join = true)
    restart_files = filter(f -> occursin("restart_t", f) && endswith(f, ".jld2"), files)

    # Helper to extract time
    function get_time_step(f)
        m = match(r"restart_t(\d+).jld2", f)
        return isnothing(m) ? 0 : parse(Int64, m[1])
    end

    sort!(restart_files, by = get_time_step)

    if length(restart_files) > keep_last_n
        files_to_delete = restart_files[1:end-keep_last_n]

        for f in files_to_delete
            # SAFETY CHECK: Do not delete if it matches the start file
            if !isempty(protect_file) && abspath(f) == abspath(protect_file)
                @info "Skipping cleanup for protected start file: $f"
                continue
            end

            rm(f; force = true)
            @info "Disk cleanup: Deleted old checkpoint $f"
        end
    end
end

end
