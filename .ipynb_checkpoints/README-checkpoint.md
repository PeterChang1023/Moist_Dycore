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
