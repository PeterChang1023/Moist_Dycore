# Optional integration check; artifacts must go into a new directory outside the repository.
# julia --project=. --compiled-modules=no --threads=1 \
#   test/validation/gcm_fixes_20261009_t42.jl /tmp/gcm_fixes_t42
using Test, JGCM, NCDatasets, TOML
include("../support/post_physics_state.jl")
const S = PostPhysicsStateSupport
root = abspath(ARGS[1])
ispath(root) && error("Validation output must be a new directory: $root")
mkpath(root)
results = Dict{String,Any}()

@testset "T42 one-day output and warm continuation" begin
    for case in (:bm_off,:bm_on)
        full_dir, split_dir = joinpath(root,"$(case)_full"),joinpath(root,"$(case)_split")
        function configuration(dir,steps)
            S.with_config(S.configuration(case; output=dir,nf=42,nlat=64,nd=20,steps);
                saving_frequency=43200,output_interval=43200,
                vars_to_output=[:u,:v,:w,:t,:q,:vor,:div,:p,:precip,:bm_dt,:bm_dq,:bm_precip],
                do_plev_output=true,pressure_levels=[85000.0,50000.0,20000.0])
        end
        full, split = configuration(full_dir,144),configuration(split_dir,72)
        full_seconds = @elapsed JGCM_Simulate(full)
        split_seconds = @elapsed JGCM_Simulate(split)
        warm = S.with_config(split; is_restart=true,
            restart_file=joinpath(split_dir,"restart/restart_t43200.jld2"))
        warm_seconds = @elapsed JGCM_Simulate(warm)
        a = Dyn_Data("full",42,43,128,64,20)
        b = Dyn_Data("split",42,43,128,64,20)
        Load_Restart_File!(a,joinpath(full_dir,"restart/restart_t86400.jld2"))
        Load_Restart_File!(b,joinpath(split_dir,"restart/restart_t86400.jld2"))
        max_error = 0.0
        for key in S.PROGNOSTICS
            @test isequal(getfield(a,key),getfield(b,key))
            max_error = max(max_error,maximum(abs,getfield(a,key)-getfield(b,key)))
        end
        for timestamp in (0,43200), suffix in ("","_plev")
            NCDataset(joinpath(full_dir,"output_t$(timestamp)$(suffix).nc"),"r") do reference
                NCDataset(joinpath(split_dir,"output_t$(timestamp)$(suffix).nc"),"r") do continued
                    for name in keys(reference)
                        @test isequal(reference[name].var[:],continued[name].var[:])
                    end
                    @test reference["sample_count"].var[:] == [72]
                    @test reference["averaging_duration"].var[:] == [43200]
                    @test reference["complete_interval"].var[:] == [1]
                    @test all(isfinite,reference["wap"].var[:])
                end
            end
        end
        results[string(case)] = Dict("steps"=>144,"timestep_seconds"=>600,
            "full_seconds"=>full_seconds,"split_cold_seconds"=>split_seconds,
            "warm_seconds"=>warm_seconds,"max_prognostic_error"=>max_error,
            "bm_active_final_step"=>any(!iszero,a.grid_bm_t_tendency),
            "min_temperature"=>minimum(a.grid_t_c),"min_pressure"=>minimum(a.grid_ps_c),
            "min_humidity"=>minimum(a.grid_q_c),"max_abs_omega"=>maximum(abs,a.grid_w_full))
    end
end
open(joinpath(root,"results.toml"),"w") do io
    TOML.print(io,results)
end
