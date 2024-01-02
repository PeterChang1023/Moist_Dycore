using JGCM
import Dates 

include("HS.jl")

#############################################################
### By CJY
start_time = Dates.now() 
###
end_day = 2
spinup_day = 0

warm_start_file_name =  "None" # "0_10day_test_warm_start_all.dat"
initial_day = 5 # warm start day

physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 65.0, "Δθ_z" => 10.0) ### 60.0
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day)

# Finalize_Output!(op_man, "5_10day_test_warm_start_final.dat", "5_10day_test_warm_start_all.dat")

# Finalize_Output!(op_man, "RH80_PR10_200day_startfrom_0day_final.dat", "RH80_PR10_200day_startfrom_0day_all.dat")

Finalize_Output!(op_man, "test_final.dat", "test_all.dat")
### time ###
final_time = Dates.now() 
all_time = final_time - start_time
@info Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(all_time))) 
############



