using Test, JGCM, JLD2
isdefined(@__MODULE__, :PostPhysicsStateSupport) || include("support/post_physics_state.jl")

@testset "Forecast restart compatibility (F4)" begin
    S = PostPhysicsStateSupport
    mktempdir() do dir
        f = S.frame(:hs; output=dir, steps=6)
        for _ in 1:3; S.step!(f); end
        contract = Restart_Contract(f.config,f.vert,f.atmo)
        @test Restart_Contract(f.config) == contract
        manager = Restart_Manager(dir,600; contract)
        Write_Restart_File(manager,f.dyn,1800)
        checkpoint = joinpath(dir,"restart_t1800.jld2")
        warm = S.with_config(f.config; is_restart=true,restart_file=checkpoint)
        restored = S.frame(warm)
        original = S.snapshot(f.dyn)
        @test all(getfield(restored.dyn,key)==value for (key,value) in original)
        @test restored.dyn.grid_w_full ≈ f.dyn.grid_w_full rtol=1e-13 atol=1e-13
        jldopen(checkpoint,"r") do file
            @test file["restart_format_version"] == 4
            @test file["restart_contract"] == contract
        end
        for changes in ((Δt=300,), (model_type=:Barotropic,), (radius=6.0e6,),
                        (omega=6e-5,), (grav=10.0,), (robert_coef=0.02,),
                        (implicit_coef=0.6,), (damping_coef=1e-5,),
                        (moisture_processes=false,), (vert_coord_option="uneven_sigma",))
            changed = S.with_config(warm; changes...)
            @test_throws ErrorException S.frame(changed)
        end
        changed_physics = merge(S.parameters(:bm_on), Dict{String,Any}("do_LRF"=>true))
        allowed = S.with_config(warm; name="new forcing", physics_params=changed_physics)
        @test Restart_Contract(allowed) == contract
        allowed_frame = S.frame(allowed)
        @test all(getfield(allowed_frame.dyn,key)==value for (key,value) in original)

        fields = filter(key -> getfield(f.dyn,key) isa Array,fieldnames(Dyn_Data))
        # Every rejection must be atomic, including invalid arrays appearing
        # after valid arrays in the file and optional diagnostics omitted safely.
        function fixture(path; version=4, metadata=contract, time=1800, omit=(), replace=Dict())
            jldopen(path,"w") do file
                file["restart_format_version"] = version
                file["saved_time"] = time
                file["dimensions"] = (f.dyn.num_fourier,f.dyn.num_spherical,
                                      f.dyn.nλ,f.dyn.nθ,f.dyn.nd)
                metadata === nothing || (file["restart_contract"] = metadata)
                for key in fields
                    key in omit && continue
                    file["state/$key"] = get(replace,key,getfield(f.dyn,key))
                end
            end
        end
        target = Dyn_Data("target",f.dyn.num_fourier,f.dyn.num_spherical,
                          f.dyn.nλ,f.dyn.nθ,f.dyn.nd)
        before = S.snapshot(target; fields)
        bad = joinpath(dir,"invalid.jld2")
        bad_ak = copy(contract["ak"]); bad_ak[2] += 1.0
        for options in ((omit=(:grid_u_p,),), (omit=(:spe_lnps_c,),),
                        (replace=Dict(:grid_q_c=>zeros(Float64,1,1,1)),),
                        (replace=Dict(:grid_t_c=>Float32.(f.dyn.grid_t_c)),),
                        (replace=Dict(:grid_u_c=>fill(NaN,size(f.dyn.grid_u_c))),),
                        (replace=Dict(:grid_ps_c=>zeros(size(f.dyn.grid_ps_c))),),
                        (replace=Dict(:grid_q_c=>fill(1.1,size(f.dyn.grid_q_c))),),
                        (metadata=merge(contract,Dict("ak"=>bad_ak)),),
                        (metadata=merge(contract,Dict("rdgas"=>300.0)),),
                        (time=-1,), (time=1800.5,), (time=1500,),
                        (time=typemax(UInt64),),
                        (metadata=nothing,))
            fixture(bad; options...)
            @test_throws ErrorException Load_Restart_File!(target,bad; expected_contract=contract)
            @test all(isequal(getfield(target,key),value) for (key,value) in before)
        end
        optional = setdiff(fields,JGCM.Restart_Manager_Module.REQUIRED_RESTART_FIELDS)
        fixture(bad; omit=optional)
        minimal = S.with_config(warm; restart_file=bad)
        minimal_frame = S.frame(minimal)
        @test all(getfield(minimal_frame.dyn,key)==value for (key,value) in original)
        @test all(isfinite,minimal_frame.dyn.grid_w_full)
        @test minimal_frame.dyn.grid_div ≈ f.dyn.grid_div

        for version in (2,3)
            fixture(bad; version,metadata=nothing)
            legacy = S.with_config(warm; restart_file=bad)
            @test_throws ErrorException S.frame(legacy)
            verified = S.with_config(legacy; restart_legacy_contract=contract)
            loaded = S.frame(verified)
            @test all(getfield(loaded.dyn,key)==value for (key,value) in original)
        end
        @test_throws ErrorException Load_Restart_File!(target,bad;
                                      expected_contract=contract,allow_partial=true)
        # Other equation sets also need reconstructed grid curl/divergence.
        two_d_config = S.with_config(f.config; model_type=:Barotropic)
        two_d_contract = Restart_Contract(two_d_config,nothing,f.atmo)
        fixture(bad; metadata=two_d_contract,omit=optional)
        two_d_warm = S.with_config(two_d_config; is_restart=true,restart_file=bad)
        Initialize_Atmos_State!(f.mesh,f.atmo,target,nothing,two_d_warm)
        @test target.grid_vor ≈ f.dyn.grid_vor
        @test target.grid_div ≈ f.dyn.grid_div
        # Continue the loaded history rather than rebuilding a cold Euler step.
        restored.integrator.time = 1800
        restored.integrator.init_step = false
        S.SI.Update_Init_Step!(restored.semi)
        for _ in 1:3
            S.step!(f); S.step!(restored)
            @test all(getfield(restored.dyn,key)==getfield(f.dyn,key) for key in S.PROGNOSTICS)
        end
    end
end
