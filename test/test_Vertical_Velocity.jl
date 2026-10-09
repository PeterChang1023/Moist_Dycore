using Test, JGCM
isdefined(@__MODULE__, :PostPhysicsStateSupport) || include("support/post_physics_state.jl")

# Independent column continuity evaluation, expressed as interface horizontal
# mass-flux divergence and the log-pressure quadrature used at full levels.
function reference_omega(f)
    d, v = f.dyn, f.vert
    gl, gt = zeros(size(d.grid_ps_c)), zeros(size(d.grid_ps_c))
    JGCM.Spectral_Dynamics_Module.Compute_Gradients!(f.mesh, d.spe_lnps_c, gl, gt)
    gl .*= d.grid_ps_c; gt .*= d.grid_ps_c
    omega = similar(d.grid_w_full)
    for j in 1:f.mesh.nθ, i in 1:f.mesh.nλ
        cumulative = 0.0
        for k in 1:f.mesh.nd
            advection = d.grid_u_c[i,j,k]*gl[i,j,1] + d.grid_v_c[i,j,k]*gt[i,j,1]
            flux_div = d.grid_div[i,j,k]*d.grid_Δp[i,j,k] + v.Δbk[k]*advection
            upper_weight = d.grid_lnp_full[i,j,k] - d.grid_lnp_half[i,j,k]
            lower_weight = d.grid_lnp_half[i,j,k+1] - d.grid_lnp_full[i,j,k]
            pressure_change = -cumulative*upper_weight - (cumulative+flux_div)*lower_weight
            pressure_advection = (v.bk[k]*upper_weight + v.bk[k+1]*lower_weight)*advection
            omega[i,j,k] = d.grid_p_full[i,j,k]*(pressure_change+pressure_advection)/d.grid_Δp[i,j,k]
            cumulative += flux_div
        end
    end
    omega
end

@testset "Accepted-state vertical pressure velocity" begin
    S = PostPhysicsStateSupport
    f = S.frame(:none)
    d = f.dyn
    fill!(d.grid_ps_c, 1e5); fill!(d.grid_lnps, log(1e5))
    Trans_Grid_To_Spherical!(f.mesh, d.grid_lnps, d.spe_lnps_c)
    S.PG.Pressure_Variables!(f.vert,d.grid_ps_c,d.grid_p_half,d.grid_Δp,
                            d.grid_lnp_half,d.grid_p_full,d.grid_lnp_full)
    fill!(d.grid_u_c,0); fill!(d.grid_v_c,0); fill!(d.grid_div,0)
    S.SD._refresh_vertical_velocity!(f.mesh,f.vert,f.atmo,d)
    @test maximum(abs,d.grid_w_full) < 1e-12
    @test all(iszero,d.grid_M_half)
    # Uniform divergence gives dp/dt = -D*p and no sigma-interface motion.
    fill!(d.grid_div,2e-6)
    before = S.snapshot(d)
    tendencies = S.snapshot(d; fields=(:grid_δu,:grid_δv,:grid_δps,:grid_δt))
    S.SD._refresh_vertical_velocity!(f.mesh,f.vert,f.atmo,d)
    @test d.grid_w_full ≈ -2e-6 .* d.grid_p_full rtol=1e-12 atol=1e-12
    @test maximum(abs,d.grid_M_half) < 1e-12
    @test all(getfield(d,key)==value for (key,value) in before)
    @test all(getfield(d,key)==value for (key,value) in tendencies)

    for case in (:none,:hs,:pbl,:bm_off,:bm_on)
        f = S.frame(case)
        for _ in 1:8
            old = reference_omega(f)
            S.step!(f)
            expected = reference_omega(f)
            @test f.dyn.grid_w_full ≈ expected rtol=1e-12 atol=1e-12
            @test maximum(abs,expected-old) > 1e-9
            @test all(iszero,f.dyn.grid_M_half[:,:,1])
            @test all(iszero,f.dyn.grid_M_half[:,:,end])
        end
    end
end
