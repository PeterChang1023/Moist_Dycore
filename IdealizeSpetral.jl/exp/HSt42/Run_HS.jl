using JGCM

include("HS.jl")

#############################################################
end_day = 2
spinup_day = 0

physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 65.0, "Δθ_z" => 10.0) ### 60.0
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day)

Finalize_Output!(op_man, "test_final.dat", "test_all.dat")



