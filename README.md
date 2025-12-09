# Moist_Dycore

## Overview

**Moist_Dycore** is an idealized atmospheric dynamical-core model with a simplified
representation of large-scale moisture and precipitation. The model is designed
for controlled process studies focusing on the interaction between atmospheric
dynamics, thermodynamics, and latent heating.

The dry backbone follows a classical hydrostatic primitive-equation dynamical
core with **Held–Suarez–type mean thermal forcing** and **Rayleigh friction**.
On top of this dry framework, a minimal large-scale condensation and
precipitation scheme is implemented to introduce moist processes while
preserving interpretability and numerical transparency.

This code is intended for **idealized numerical experiments and educational
use**, not for operational weather prediction or comprehensive climate modeling.

---

## Purpose and Scope

The main objectives of this project are:

- To provide a lightweight and readable atmospheric dynamical core  
- To study the impact of moisture and latent heat release on large-scale
  circulation  
- To serve as a testbed for numerical methods and simplified physical
  parameterizations  

Typical users include atmospheric science students, researchers conducting
idealized simulations, and developers experimenting with numerical or physical
model components.

---

## Key Features

- Global hydrostatic primitive-equation dynamical core  
- Held–Suarez-style Newtonian temperature relaxation  
- Rayleigh damping near the lower boundary  
- Prognostic specific humidity  
- Large-scale condensation and precipitation  
- Latent heat release coupled to the thermodynamic equation  
- Modular and transparent code design  

---

## Repository Structure
- README.md
- IdealizeSpetral
    - exp
        - HSt42 (Run Dycore at here)
        - HSt21
        - Barotropic
        - Shallow_Water
    - src
        - Atmos_Spectral
            - Spectral_Dynamics.jl (this is the main code of Dycore)
            - Spectral_Spherical_Mesh.jl
            - Output_Manager.jl
            - Time_Integrator.jl
            - Vert_Coordinate.jl
            - Shallow_Water_Dynamics.jl
            - Semi_Implicit.jl
            - Press_And_Geopot.jl
            - Gauss_And_Legendre.jl
            - Dyn_Data.jl
            - Barotropic_Dynamics.jl
            - Atmo_Data.jl
        - Atmos_Param
            - HS_Forcing.jl
    - julia
    - Manifest.toml
    - Project.toml

---
## User Guide: Running a Held–Suarez Test Case (HSt42)

This guide shows the minimum steps required to set up the environment and
run a test Held–Suarez (HS) experiment.  
Please follow the steps **in order**.

---

### Step 1. Create the Conda Environment (once)

```
conda env create -f Dycore_environment.yml
```
This creates a Conda environment named Dycore_env.

### Step 2. Activate the Environment (important)

```
conda activate Dycore_env
```
⚠️ Important:
All subsequent commands (Julia, JupyterLab, Git) must be executed after
activating Dycore_env.

### Step 3. enter main fold
```
cd IdealizeSpetral
```

### Step 4: find the location of Python, which is useful for 
```
which python
---> output: ~miniconda3/envs/Dycore_env/bin/python
(copy output, it will be used at Step 6.)
```

### Step 5: call julia
```
julia
```

### Step 6: in Julia, using ... (similar to import in Python)
```
julia> using Pkg
```

### Step 7: point the python location for Julia 
```
julia> ENV["PYTHON"] = "~miniconda3/envs/Dycore_env/bin/python"
```

### Step 8: import PyCall
```
julia> Pkg.build("PyCall")
```

### Step 9: import PyPlot
```
julia> Pkg.build("PyPlot")
```

### Step 10: enter Julia package mode
```
julia> ]
```

### Step 11: install all packages
```
pkg> dev . 
```

### Step 12: precompile the all packages in JGCM
```
(backspace back to julia>)
julia> using JGCM 
```

### Step 13: exit Julia and run the test run
```
ctrl + D
cd exp/HSt42
julia For_test_Run_HS.jl
```

### Notion edition:
https://www.notion.so/29c9c6fdeb6c80b0855aca2c848ff751?v=29c9c6fdeb6c803889e4000cbb063880&p=2c49c6fdeb6c80ae9239c782d351d6ba&pm=s



