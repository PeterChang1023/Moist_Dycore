using Test, JGCM, NCDatasets
isdefined(@__MODULE__, :PostPhysicsStateSupport) || include("support/post_physics_state.jl")

@testset "Averaging windows and restart phase (F9)" begin
    S = PostPhysicsStateSupport
    OM = JGCM.Output_Manager_Module
    f = S.frame(:hs)
    function manager(dir; start=0,spinup=0.0)
        Output_Manager(f.mesh,f.vert,f.atmo,start,start+10800,[:u];
            filename=joinpath(dir,"output.nc"),output_interval=3600,
            timestep=600,spinup_day=spinup)
    end
    function sample!(mgr,time)
        fill!(f.dyn.grid_u_c,time)
        Update_Output!(mgr,f.dyn,time)
    end
    function records(path)
        NCDataset(path,"r") do ds
            @test ds["time"].attrib["bounds"] == "time_bounds"
            @test ds["ua"].attrib["cell_methods"] == "time: mean"
            (time=ds["time"].var[:],bounds=ds["time_bounds"].var[:,:],
             count=ds["sample_count"].var[:],duration=ds["averaging_duration"].var[:],
             complete=ds["complete_interval"].var[:],u=ds["ua"].var[1,1,1,:])
        end
    end
    mktempdir() do dir
        mgr = manager(dir)
        for time in 600:600:5400; sample!(mgr,time); end
        before = copy(mgr.acc_raw[:u])
        @test_throws ArgumentError OM.Rotate_NC_Chunk!(mgr,5400)
        @test mgr.acc_raw[:u] == before
        @test mgr.sample_counter == 3
        @test !isfile(joinpath(dir,"output_t5400.nc"))
        for time in 6000:600:7200; sample!(mgr,time); end
        OM.Rotate_NC_Chunk!(mgr,7200)
        for time in 7800:600:10800; sample!(mgr,time); end
        Finalize_Output!(mgr)
        a,b = records(joinpath(dir,"output_t0.nc")),records(joinpath(dir,"output_t7200.nc"))
        @test vcat(a.time,b.time) ≈ [3600,7200,10800]./86400
        @test hcat(a.bounds,b.bounds) ≈ [0 3600 7200;3600 7200 10800]./86400
        @test vcat(a.count,b.count) == [6,6,6]
        @test vcat(a.duration,b.duration) == [3600,3600,3600]
        @test vcat(a.complete,b.complete) == [1,1,1]
        @test vcat(a.u,b.u) == [2100,5700,9300]
    end
    mktempdir() do dir
        mgr = manager(dir)
        for time in 600:600:5400; sample!(mgr,time); end
        Finalize_Output!(mgr)
        r = records(joinpath(dir,"output_t0.nc"))
        @test r.bounds ≈ [0 3600;3600 5400]./86400
        @test r.count == [6,3]
        @test r.duration == [3600,1800]
        @test r.complete == [1,0]
        @test r.u == [2100,4800]
    end
    mktempdir() do dir
        mgr = manager(dir)
        for time in 600:600:3600; sample!(mgr,time); end
        Finalize_Output!(mgr)
        warm = manager(dir;start=3600)
        for time in 4200:600:7200; sample!(warm,time); end
        Finalize_Output!(warm)
        a,b = records(joinpath(dir,"output_t0.nc")),records(joinpath(dir,"output_t3600.nc"))
        @test vcat(a.u,b.u) == [2100,5700]
        @test hcat(a.bounds,b.bounds) ≈ [0 3600;3600 7200]./86400
        @test vcat(a.count,b.count) == [6,6]
    end
    mktempdir() do dir
        @test_throws ArgumentError manager(dir;start=1800)
        @test_throws ArgumentError manager(dir;spinup=1800/86400)
        mgr = manager(dir;spinup=3600/86400)
        for time in 600:600:7200; sample!(mgr,time); end
        @test_throws ArgumentError Update_Output!(mgr,f.dyn,8400)
        @test_throws ArgumentError Update_Output!(mgr,f.dyn,7200)
        Finalize_Output!(mgr)
        r = records(joinpath(dir,"output_t0.nc"))
        @test r.bounds ≈ reshape([3600,7200]./86400,2,1)
        @test r.count == [6]
        @test r.complete == [1]
    end
    for changed in ((saving_frequency=5400,), (spinup_day=1800/86400,))
        cfg = S.with_config(S.configuration(:hs);output_interval=3600,changed...)
        @test_throws ArgumentError JGCM.Driver.validate_config(cfg)
    end
    mktempdir() do dir
        source = S.frame(:hs)
        for _ in 1:3; S.step!(source); end
        Write_Restart_File(Restart_Manager(dir,0;contract=Restart_Contract(source.config,source.vert,source.atmo)),
                           source.dyn,1800)
        run = joinpath(dir,"warm")
        cfg = S.with_config(S.configuration(:hs;output=run,steps=1);
            output_interval=3600,is_restart=true,restart_file=joinpath(dir,"restart_t1800.jld2"))
        @test_throws ArgumentError JGCM_Simulate(cfg)
        @test !isfile(joinpath(run,"output_t1800.nc"))
    end
end
