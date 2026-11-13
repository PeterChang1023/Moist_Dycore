module Output_Manager_Module

using Base
using Dates
using NCDatasets

using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Vert_Coordinate_Module
using ..Vertical_Interpolation_Module
using ..Output_Mappings_Module
using ..Variable_Mappings_Module

export Output_Manager, Update_Output!, Finalize_Output!, Rotate_NC_Chunk!
export PrimitiveEquationMode, BarotropicMode, ShallowWaterMode

# ==============================================================================
# 1. Type Hierarchy (The Backbone of Dispatch)
# ==============================================================================
abstract type AbstractModelMode end

struct PrimitiveEquationMode <: AbstractModelMode end
struct BarotropicMode <: AbstractModelMode end
struct ShallowWaterMode <: AbstractModelMode end

# Helper to map legacy symbols to Types
function Mode_Factory(sym::Symbol)
    if sym == :PrimitiveEquation
        return PrimitiveEquationMode()
    elseif sym == :Barotropic
        return BarotropicMode()
    elseif sym == :ShallowWater
        return ShallowWaterMode()
    else
        error("Output_Manager: Unknown model_mode symbol :$sym")
    end
end

# Helper to map Types back to Symbols (for file attributes/legacy APIs)
mode_symbol(::PrimitiveEquationMode) = :PrimitiveEquation
mode_symbol(::BarotropicMode)        = :Barotropic
mode_symbol(::ShallowWaterMode)      = :ShallowWater

# Only inspect this output series; unrelated NetCDF products can share a directory.
function _output_series_files(filename::AbstractString)
    base = splitext(filename)[1]
    dir, stem = dirname(base), basename(base)
    isdir(dir) || return Pair{String,Union{Nothing,Int64}}[]
    matches = Pair{String,Union{Nothing,Int64}}[]
    for name in readdir(dir)
        if name in (stem * ".nc", stem * "_plev.nc")
            push!(matches, joinpath(dir, name) => nothing)
        elseif startswith(name, stem * "_t") && endswith(name, ".nc")
            suffix = chop(name; head=length(stem) + 2, tail=3)
            endswith(suffix, "_plev") && (suffix = chop(suffix; tail=5))
            occursin(r"^\d+$", suffix) || continue
            timestamp = tryparse(Int64, suffix)
            timestamp === nothing || push!(matches, joinpath(dir, name) => timestamp)
        end
    end
    return matches
end

function validate_output_series(filename::AbstractString, start_time::Integer;
                                day_to_sec::Integer=86400)
    for (path, chunk_start) in _output_series_files(filename)
        if start_time == 0 || chunk_start === nothing || chunk_start >= start_time
            throw(ArgumentError("existing output conflicts with this run: $path; " *
                                "use a fresh output filename/directory or archive the old series"))
        end
        # A chunk whose filename predates the restart may still contain later records.
        NCDataset(path, "r") do ds
            haskey(ds, "time") || throw(ArgumentError("output has no time coordinate: $path"))
            times = ds["time"].var[:]
            get(ds["time"].attrib, "units", "") == "days since 2000-01-01 00:00:00" ||
                throw(ArgumentError("cannot establish output time coverage: $path"))
            all(t -> isfinite(t) && t <= start_time / day_to_sec, times) ||
                throw(ArgumentError("output extends beyond restart time: $path; use a fresh output directory"))
        end
    end
    return nothing
end

function _check_chunk_paths(base, timestamp)
    for path in ("$(base)_t$(timestamp).nc", "$(base)_t$(timestamp)_plev.nc")
        ispath(path) && throw(ArgumentError("refusing to overwrite existing output: $path"))
    end
end

# ==============================================================================
# 2. The Parameterized Struct
# ==============================================================================
mutable struct Output_Manager{M <: AbstractModelMode}
    # Grid Metadata
    nλ::Int64
    nθ::Int64
    nd::Int64

    # --- Physics Context ---
    atmo_data::Atmo_Data
    mode::M

    # --- Configuration ---
    do_plev_output::Bool             # true = also write pressure-level interpolated output
    target_levels::Vector{Float64}
    log_targets::Vector{Float64}

    # --- Buffers (Pre-allocated, 3D only) ---
    interpolation_cache::Pressure_Interpolation_Cache
    interp_buffer::Array{Float64, 3}

    vert_ak::Vector{Float64}
    vert_bk::Vector{Float64}

    # --- Time Control ---
    day_to_sec::Int64
    start_time::Int64
    current_time::Int64
    spinup_time::Int64
    output_interval::Int64
    sample_counter::Int64
    output_index::Int64
    window_start_time::Int64
    sampling_interval::Int64

    # --- File Handles ---
    ds_plev::Union{NCDataset, Nothing}  # pressure-level output (optional, 3D only)
    ds_raw::Union{NCDataset, Nothing}   # native sigma/hybrid (or 2D native) — always created

    # --- Accumulators ---
    # The native-grid mean is shared by raw and pressure-level output. Keeping
    # one accumulator avoids adding every model field twice on every time step.
    acc_raw::Dict{Symbol, Array{Float64}}

    active_symbols::Vector{Symbol}
    var_info_map::Dict{Symbol, VarMeta}

    # --- Chunk I/O State ---
    base_filename::String            # base path without extension, e.g. "exp/HSt42/output"
    mesh_λc::Vector{Float64}        # longitude grid in radians (for file re-creation on rotation)
    mesh_θc::Vector{Float64}        # latitude grid in radians
    vert_coord_cache::Any           # Vert_Coordinate or nothing
    global_meta::Dict{String, Any}  # NC global attributes
end

# ==============================================================================
# 3. Internal Initialization Logic
# ==============================================================================
function _Init_Single_File(
    filename, mesh, var_info_map,
    requested_vars, file_type::Symbol, target_levels=Float64[];
    vert_coord=nothing,
    global_meta::Dict{String, Any}=Dict{String, Any}()
)
    ispath(filename) && throw(ArgumentError("refusing to overwrite existing output: $filename"))
    ds = NCDataset(filename, "c")

    # --- Global Attributes (CF-1.11 Compliance) --- #
    for (key, val) in global_meta
        ds.attrib[key] = val
    end
    ds.attrib["Conventions"] = "CF-1.11"
    ds.attrib["history"]     = "Created $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))"
    ds.attrib["source"]      = "idealized-spectral-gcm"
    # --- Global Attributes (CF-1.11 Compliance) --- #

    # --- Dimensions --- #
    defDim(ds, "lon",  mesh.nλ)
    defDim(ds, "lat",  mesh.nθ)
    defDim(ds, "time", Inf)
    # --- Dimensions --- #

    # --- Coordinate Variables & Vertical Metadata --- #
    # Longitude
    v_lon = defVar(ds, "lon", Float64, ("lon",),
                   attrib=Dict("units"=>"degrees_east", "standard_name"=>"longitude"))
    v_lon[:] = mesh.λc * 180/pi

    # Latitude
    v_lat = defVar(ds, "lat", Float64, ("lat",),
                   attrib=Dict("units"=>"degrees_north", "standard_name"=>"latitude"))
    v_lat[:] = mesh.θc * 180/pi

    # Time
    defVar(ds, "time", Float64, ("time",),
           attrib=Dict("units"=>"days since 2000-01-01 00:00:00", "calendar"=>"proleptic_gregorian",
                       "bounds"=>"time_bounds"))
    defDim(ds, "bounds", 2)
    defVar(ds, "time_bounds", Float64, ("bounds", "time"), attrib=Dict(
        "units"=>"days since 2000-01-01 00:00:00", "calendar"=>"proleptic_gregorian"))
    defVar(ds, "sample_count", Int64, ("time",), attrib=Dict(
        "units"=>"1", "long_name"=>"number of equally spaced endpoint samples in time mean"))
    defVar(ds, "averaging_duration", Float64, ("time",), attrib=Dict(
        "units"=>"s", "long_name"=>"duration of averaging window"))
    defVar(ds, "complete_interval", Int8, ("time",), attrib=Dict(
        "long_name"=>"whether the record covers a complete requested output interval",
        "flag_values"=>Int8[0,1], "flag_meanings"=>"partial complete"))

    if file_type == :plev
        defDim(ds, "plev", length(target_levels))
        v_lev = defVar(ds, "plev", Float64, ("plev",),
                       attrib=Dict("units"=>"Pa", "standard_name"=>"air_pressure", "positive"=>"down"))
        v_lev[:] = target_levels

    elseif file_type == :raw && vert_coord !== nothing
        # Isca-style native vertical metadata. These one-dimensional pressure
        # axes are approximate reference values; exact column pressure is
        # reconstructed from interface pk, bk, and the output-record ps.
        defDim(ds, "pfull", mesh.nd)
        defDim(ds, "phalf", mesh.nd + 1)

        reference_ps = fill(Float64(vert_coord.p_ref), 1, 1)
        reference_pfull = zeros(Float64, 1, 1, mesh.nd)
        Compute_Pressure_Grid!(reference_pfull, vert_coord.ak, vert_coord.bk, reference_ps)
        reference_phalf = vert_coord.ak .+ vert_coord.bk .* vert_coord.p_ref

        v_pfull = defVar(ds, "pfull", Float64, ("pfull",), attrib=Dict(
            "units"=>"hPa", "long_name"=>"approx full pressure level",
            "positive"=>"down", "axis"=>"Z"))
        v_pfull[:] = vec(reference_pfull) .* 0.01

        v_phalf = defVar(ds, "phalf", Float64, ("phalf",), attrib=Dict(
            "units"=>"hPa", "long_name"=>"approx half pressure level",
            "positive"=>"down", "axis"=>"Z"))
        v_phalf[:] = reference_phalf .* 0.01

        defVar(ds, "pk", Float64, ("phalf",), attrib=Dict(
            "units"=>"Pa", "long_name"=>"vertical coordinate pressure values"))[:] = vert_coord.ak
        defVar(ds, "bk", Float64, ("phalf",), attrib=Dict(
            "units"=>"1", "long_name"=>"vertical coordinate sigma values"))[:] = vert_coord.bk
    end
    # --- Coordinate Variables & Vertical Metadata --- #

    # --- Data Variables --- #
    for sym in requested_vars
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]

        dims = if meta.dims == 3
            (file_type == :plev) ? ("lon", "lat", "plev", "time") : ("lon", "lat", "pfull", "time")
        else
            ("lon", "lat", "time")
        end

        v = defVar(ds, meta.nc_name, Float64, dims)
        v.attrib["units"]         = meta.units
        v.attrib["long_name"]     = meta.long_name
        v.attrib["standard_name"] = meta.std_name
        v.attrib["cell_methods"]  = "time: mean"
    end
    # --- Data Variables --- #

    return ds
end

# ==============================================================================
# 4. The Constructor
# ==============================================================================
function Output_Manager(
    mesh, vert_coord, atmo_data,
    start_time, ::Any, requested_vars;
    filename::String = "output.nc",
    do_plev_output::Bool = false,
    pressure_levels::Vector{Float64} = Float64[],
    day_to_sec::Int64 = 86400,
    output_interval::Int64 = 86400,
    timestep::Int64 = 0,
    spinup_day::Float64 = 0.0,
    model_mode::Symbol = :PrimitiveEquation,
    institute::String = "My Research Lab",
    experiment_id::String = "Simulation_v1"
)
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    output_interval > 0 || throw(ArgumentError("output_interval must be positive"))
    spinup_day >= 0 || throw(ArgumentError("spinup_day must be non-negative"))
    start_time >= 0 && start_time % output_interval == 0 || throw(
        ArgumentError("output start/restart time must align with output_interval; no partial accumulator is restored"))
    spinup_seconds = spinup_day * day_to_sec
    isfinite(spinup_seconds) && isinteger(spinup_seconds) || throw(
        ArgumentError("spinup duration must be a finite integer number of seconds"))
    Int64(spinup_seconds) % output_interval == 0 || throw(
        ArgumentError("spinup duration must align with output_interval"))
    timestep >= 0 && (timestep == 0 || output_interval % timestep == 0) || throw(
        ArgumentError("timestep must be non-negative and divide output_interval"))
    validate_output_series(filename, start_time; day_to_sec)
    if do_plev_output
        !isempty(pressure_levels) || throw(
            ArgumentError("pressure-level output requires at least one pressure level"),
        )
        all(p -> isfinite(p) && p > 0, pressure_levels) || throw(
            ArgumentError("pressure levels must be positive and finite"),
        )
    end

    # --- Instantiate the Mode Type --- #
    mode_obj    = Mode_Factory(model_mode)
    if do_plev_output && !isa(mode_obj, PrimitiveEquationMode)
        throw(ArgumentError("pressure-level output is only available for PrimitiveEquation"))
    end
    nλ, nθ, nd  = mesh.nλ, mesh.nθ, mesh.nd
    spinup_time = start_time + Int64(spinup_day * day_to_sec)
    active_symbols = unique(copy(requested_vars))

    # Derive base filename (strip ".nc" extension)
    base_fn = splitext(filename)[1]

    # Global Metadata Dictionary for the NetCDF header
    global_meta = Dict(
        "title"       => "Output of idealized-spectral-gcm",
        "institution" => institute,
        "experiment"  => experiment_id,
        "references"  => "https://github.com/Brownian-Motion-99/idealized-spectral-gcm.git",
        "gravity" => atmo_data.grav,
        "gravity_units" => "m s-2",
        "dry_air_gas_constant" => atmo_data.rdgas,
        "dry_air_gas_constant_units" => "J kg-1 K-1",
        "output_interval_seconds" => output_interval,
        "averaging_sample_convention" => "equally spaced model-step endpoints; time bounds enclose preceding steps",
    )

    # Define Buffers (Initialize as empty by default)
    interpolation_cache = Pressure_Interpolation_Cache(0, 0, 1, 0)
    interp_buf = zeros(Float64, 0, 0, 0)
    vert_ak    = Float64[]
    vert_bk    = Float64[]
    # --- Instantiate the Mode Type --- #

    # --- 3D Specific Initialization --- #
    if isa(mode_obj, PrimitiveEquationMode)
        # Validation
        if !(:ps in active_symbols); push!(active_symbols, :ps); end
        if (:z in active_symbols) && !(:t in active_symbols); error("Need :t for :z"); end

        if vert_coord === nothing
             error("Output_Manager: vert_coord is required for :PrimitiveEquation mode")
        end

        # Allocation (Only happens for 3D)
        n_plev     = length(pressure_levels)
        if do_plev_output
            interpolation_cache = Pressure_Interpolation_Cache(nλ, nθ, nd, n_plev)
            interp_buf = zeros(Float64, nλ, nθ, n_plev)
        end

        vert_ak = copy(vert_coord.ak)
        vert_bk = copy(vert_coord.bk)
    end

    # Fetch variable info
    var_info_map = Get_Var_Info(Val(mode_symbol(mode_obj)))
    unknown_symbols = setdiff(active_symbols, keys(var_info_map))
    isempty(unknown_symbols) || throw(
        ArgumentError("unsupported output variables for $(mode_symbol(mode_obj)): $unknown_symbols"),
    )
    # --- 3D Specific Initialization --- #

    # --- Initialize Accumulators --- #
    # acc_raw: always allocated (native grid accumulator for all model types)
    acc_raw = Dict{Symbol, Array{Float64}}()
    for sym in active_symbols
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]
        acc_raw[sym] = (meta.dims == 3) ? zeros(Float64, nλ, nθ, nd) : zeros(Float64, nλ, nθ)
    end

    # --- Initialize Accumulators --- #

    # --- Initialize Files --- #
    # ds_raw: always created for all model types (primary/native output)
    raw_file_type = isa(mode_obj, PrimitiveEquationMode) ? :raw : :simple_2d
    ds_raw = _Init_Single_File(
        "$(base_fn)_t$(start_time).nc", mesh,
        var_info_map, active_symbols, raw_file_type;
        vert_coord=vert_coord, global_meta=global_meta
    )

    # ds_plev: only for 3D with do_plev_output = true
    ds_plev = nothing
    if do_plev_output && isa(mode_obj, PrimitiveEquationMode)
        ds_plev = _Init_Single_File(
            "$(base_fn)_t$(start_time)_plev.nc", mesh,
            var_info_map, active_symbols, :plev, pressure_levels;
            vert_coord=vert_coord, global_meta=global_meta
        )
    end
    # --- Initialize Files --- #

    return Output_Manager(
        nλ, nθ, nd, atmo_data, mode_obj,
        do_plev_output, pressure_levels, log.(max.(pressure_levels, 1.0)),
        interpolation_cache, interp_buf, vert_ak, vert_bk,
        day_to_sec, start_time, start_time, spinup_time, output_interval, 0, 1,
        spinup_time, timestep,
        ds_plev, ds_raw, acc_raw, active_symbols, var_info_map,
        base_fn, collect(mesh.λc), collect(mesh.θc), vert_coord, global_meta
    )
end

# ==============================================================================
# 5. Runtime Logic (Dispatch Enabled)
# ==============================================================================

# Generic Wrapper
function Update_Output!(manager::Output_Manager{M}, dyn_data, current_time::Int64) where M
    if current_time == manager.current_time && current_time == manager.start_time
        return
    end
    elapsed = current_time - manager.current_time
    elapsed > 0 || throw(ArgumentError("output sample times must increase strictly"))
    if manager.sampling_interval == 0
        manager.output_interval % elapsed == 0 || throw(
            ArgumentError("sampling interval must divide output_interval"))
        manager.sampling_interval = elapsed
    end
    elapsed == manager.sampling_interval || throw(
        ArgumentError("output samples must be equally spaced; expected $(manager.sampling_interval) s, got $elapsed s"))
    manager.current_time = current_time
    if current_time <= manager.spinup_time; return; end

    # Fetch Data using the Mode Type (cleaner dispatch)
    live_data = Get_Dyn_Var_Map(dyn_data, Val(mode_symbol(manager.mode)))

    # CALL SPECIALIZED ACCUMULATION
    _accumulate_core!(manager.mode, manager, live_data)

    manager.sample_counter += 1

    time_elapsed = current_time - manager.start_time
    if time_elapsed > 0 && (time_elapsed % manager.output_interval == 0)
        Flush_to_Disk!(manager)
    end
end

# --- Specialized Accumulation: Primitive Equation ---
function _accumulate_core!(::PrimitiveEquationMode, mgr, live_data)
    for sym in mgr.active_symbols
        if !haskey(live_data, sym); continue; end

        raw_data = live_data[sym]

        # acc_raw: always accumulate (native sigma/hybrid grid → ds_raw)
        if haskey(mgr.acc_raw, sym)
            _accumulate_buffer!(mgr.acc_raw[sym], raw_data)
        end
    end
end

# --- Specialized Accumulation: 2D Modes ---
function _accumulate_core!(::AbstractModelMode, mgr, live_data)
    # Generic fallback for Barotropic/ShallowWater: use acc_raw
    for sym in mgr.active_symbols
        if !haskey(live_data, sym); continue; end

        if haskey(mgr.acc_raw, sym)
            mgr.acc_raw[sym] .+= live_data[sym]
        end
    end
end

# Helper to handle singleton dimensions (e.g. Surface Pressure in 3D array)
@inline function _accumulate_buffer!(buffer, data)
    if ndims(buffer) == ndims(data)
        buffer .+= data
    elseif ndims(buffer) == 2 && ndims(data) == 3
        buffer .+= view(data, :, :, 1)
    end
end

# ==============================================================================
# 6. IO Logic (Dispatch Enabled)
# ==============================================================================

function Flush_to_Disk!(manager::Output_Manager{M}; allow_partial::Bool=false) where M
    if manager.sample_counter == 0; return; end

    t_val = manager.current_time / manager.day_to_sec
    duration = manager.current_time - manager.window_start_time
    duration == manager.sample_counter * manager.sampling_interval ||
        error("Output samples do not cover the recorded averaging window")
    complete = duration == manager.output_interval
    (complete || (allow_partial && 0 < duration < manager.output_interval)) ||
        throw(ArgumentError("cannot flush an incomplete output interval except at finalization"))

    # Write time index to whichever files are open
    for ds in (manager.ds_raw,manager.ds_plev)
        ds === nothing && continue
        idx = manager.output_index
        ds["time"][idx] = t_val
        ds["time_bounds"][:,idx] = [manager.window_start_time / manager.day_to_sec,t_val]
        ds["sample_count"][idx] = manager.sample_counter
        ds["averaging_duration"][idx] = duration
        ds["complete_interval"][idx] = Int8(complete)
    end

    # Dispatch to specific writing logic
    _write_core!(manager.mode, manager)

    manager.sample_counter = 0
    manager.window_start_time = manager.current_time
    manager.output_index += 1
end

# --- Writer: Primitive Equation ---
function _write_core!(::PrimitiveEquationMode, mgr)
    t_idx        = mgr.output_index
    sample_count = Float64(mgr.sample_counter)
    var_info_map = mgr.var_info_map

    # Finalize the native-grid mean once. Both output products consume this
    # same record, matching Isca's average-then-interpolate ordering.
    for accumulator in values(mgr.acc_raw)
        accumulator ./= sample_count
    end

    # =========================================================================
    # PRIMARY: Native sigma/hybrid output (acc_raw → ds_raw, always)
    # =========================================================================
    for sym in keys(mgr.acc_raw)
        haskey(var_info_map, sym) || continue
        meta     = var_info_map[sym]
        nc_var   = mgr.ds_raw[meta.nc_name]
        raw_mean = mgr.acc_raw[sym]

        if ndims(nc_var) == 4  # (lon, lat, pfull, time)
            nc_var[:, :, :, t_idx] = raw_mean
        elseif ndims(nc_var) == 3  # (lon, lat, time)
            nc_var[:, :, t_idx] = (ndims(raw_mean) == 3) ? view(raw_mean, :, :, 1) : raw_mean
        end
    end
    NCDatasets.sync(mgr.ds_raw)

    # =========================================================================
    # OPTIONAL: Pressure-level interpolated output (same mean → ds_plev)
    # =========================================================================
    if mgr.do_plev_output && mgr.ds_plev !== nothing
        ps_avg_2d = mgr.acc_raw[:ps]
        Prepare_Interpolation!(
            mgr.interpolation_cache,
            mgr.vert_ak,
            mgr.vert_bk,
            ps_avg_2d,
            mgr.log_targets,
        )
        t_ref = get(mgr.acc_raw, :t, nothing)

        for sym in keys(mgr.acc_raw)
            haskey(var_info_map, sym) || continue
            meta        = var_info_map[sym]
            nc_var      = mgr.ds_plev[meta.nc_name]
            native_mean = mgr.acc_raw[sym]

            if ndims(nc_var) == 4  # (lon, lat, plev, time)
                Apply_Interpolation!(
                    mgr.interp_buffer,
                    native_mean,
                    mgr.interpolation_cache,
                    sym,
                    mgr.atmo_data,
                    t_ref,
                )
                nc_var[:, :, :, t_idx] = mgr.interp_buffer
            else
                nc_var[:, :, t_idx] = (ndims(native_mean) == 3) ? view(native_mean, :, :, 1) : native_mean
            end
        end
        NCDatasets.sync(mgr.ds_plev)
    end

    for accumulator in values(mgr.acc_raw)
        fill!(accumulator, 0.0)
    end
end

# --- Writer: 2D Modes (Barotropic/ShallowWater) ---
function _write_core!(::AbstractModelMode, mgr)
    t_idx        = mgr.output_index
    sample_count = Float64(mgr.sample_counter)
    var_info_map = mgr.var_info_map

    for sym in keys(mgr.acc_raw)
        if !haskey(var_info_map, sym); continue; end
        nc_name = var_info_map[sym].nc_name
        accumulator = mgr.acc_raw[sym]
        accumulator ./= sample_count
        mgr.ds_raw[nc_name][:, :, t_idx] = accumulator
        fill!(accumulator, 0.0)
    end
    NCDatasets.sync(mgr.ds_raw)
end

# ==============================================================================
# 7. Chunk Rotation
# ==============================================================================

"""
    _open_chunk_files!(manager, chunk_time)

Internal helper: open a new pair of NC output files for the given chunk start time.
Uses a NamedTuple mesh proxy so the full mesh object need not be stored in the manager.
"""
function _open_chunk_files!(manager::Output_Manager, chunk_time::Int64)
    var_info_map = manager.var_info_map

    # NamedTuple acts as a mesh proxy (supports dot-notation like a struct)
    mesh_proxy = (
        nλ = manager.nλ,
        nθ = manager.nθ,
        nd = manager.nd,
        λc = manager.mesh_λc,
        θc = manager.mesh_θc
    )

    # ds_raw: always (primary output for all model types)
    raw_file_type = isa(manager.mode, PrimitiveEquationMode) ? :raw : :simple_2d
    manager.ds_raw = _Init_Single_File(
        "$(manager.base_filename)_t$(chunk_time).nc", mesh_proxy,
        var_info_map, manager.active_symbols, raw_file_type;
        vert_coord=manager.vert_coord_cache, global_meta=manager.global_meta
    )

    # ds_plev: optional (3D with do_plev_output = true only)
    if manager.do_plev_output && isa(manager.mode, PrimitiveEquationMode)
        manager.ds_plev = _Init_Single_File(
            "$(manager.base_filename)_t$(chunk_time)_plev.nc", mesh_proxy,
            var_info_map, manager.active_symbols, :plev, manager.target_levels;
            vert_coord=manager.vert_coord_cache, global_meta=manager.global_meta
        )
    end
end

"""
    Rotate_NC_Chunk!(manager, new_chunk_time)

Close the current NC output files at a completed averaging boundary, and open fresh
time-stamped files for the next chunk. Called at each saving_frequency boundary
in the driver, coordinated with JLD2 checkpoint writes.
"""
function Rotate_NC_Chunk!(manager::Output_Manager, new_chunk_time::Int64)
    _check_chunk_paths(manager.base_filename, new_chunk_time)
    new_chunk_time == manager.current_time && new_chunk_time % manager.output_interval == 0 ||
        throw(ArgumentError("chunk rotation must occur at the current completed output boundary"))
    manager.sample_counter == 0 || throw(
        ArgumentError("cannot rotate with a partial averaging window"))

    # Close current chunk files
    if manager.ds_plev !== nothing
        close(manager.ds_plev)
        manager.ds_plev = nothing
    end
    if manager.ds_raw !== nothing
        close(manager.ds_raw)
        manager.ds_raw = nothing
    end

    # Open next chunk files
    _open_chunk_files!(manager, new_chunk_time)

    # Reset time index for the new file
    manager.output_index = 1

    @info "NC output rotated: new chunk at t=$(new_chunk_time)s → $(manager.base_filename)_t$(new_chunk_time).nc"
end

# ==============================================================================
# 8. Finalization
# ==============================================================================

function Finalize_Output!(manager::Output_Manager)
    if manager.sample_counter > 0
        Flush_to_Disk!(manager; allow_partial=true)
    end
    if manager.ds_plev !== nothing; close(manager.ds_plev); end
    if manager.ds_raw  !== nothing; close(manager.ds_raw);  end
end

end
