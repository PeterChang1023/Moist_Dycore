using Test, JGCM, NCDatasets
isdefined(@__MODULE__, :PostPhysicsStateSupport) || include("support/post_physics_state.jl")

@testset "Output series integrity (F3)" begin
    OM = JGCM.Output_Manager_Module
    S = PostPhysicsStateSupport
    for conflict in ("output_t0.nc", "output_t999999.nc", "output_t999999_plev.nc",
                     "output.nc", "restart/restart_t600.jld2")
        mktempdir() do dir
            path = joinpath(dir, conflict)
            mkpath(dirname(path))
            write(path, "old run bytes")
            write(joinpath(dir, "progress.log"), "old log")
            cfg = S.configuration(:hs; output=dir, steps=1)
            @test_throws ArgumentError JGCM_Simulate(cfg)
            @test read(path, String) == "old run bytes"
            @test read(cfg.logger, String) == "old log"
        end
    end
    mktempdir() do dir
        write(joinpath(dir, "unrelated.nc"), "keep")
        f = S.frame(:hs)
        filename = joinpath(dir, "output.nc")
        manager = Output_Manager(f.mesh, f.vert, f.atmo, 0, 1200, [:u];
                                 filename, output_interval=600)
        Update_Output!(manager, f.dyn, 600)
        # Both native and pressure-level target collisions must be caught before flush/close.
        conflict = joinpath(dir, "output_t600_plev.nc")
        write(conflict, "old pressure output")
        @test_throws ArgumentError OM.Rotate_NC_Chunk!(manager, 600)
        @test read(conflict, String) == "old pressure output"
        @test manager.ds_raw !== nothing
        Finalize_Output!(manager)
        rm(conflict)
        @test_throws ArgumentError OM.validate_output_series(filename, 300)
        @test OM.validate_output_series(filename, 600) === nothing
        warm = Output_Manager(f.mesh, f.vert, f.atmo, 600, 1200, [:u];
                              filename, output_interval=600)
        Update_Output!(warm, f.dyn, 1200)
        Finalize_Output!(warm)
        @test_throws ArgumentError OM.validate_output_series(filename, 600)
        @test read(joinpath(dir, "unrelated.nc"), String) == "keep"
        before = read(joinpath(dir, "output_t0.nc"))
        @test_throws ArgumentError OM._Init_Single_File(joinpath(dir, "output_t0.nc"),
            f.mesh, Dict(), Symbol[], :raw)
        @test read(joinpath(dir, "output_t0.nc")) == before
    end
end
