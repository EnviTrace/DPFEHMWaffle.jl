# DPFEHMWaffle.jl
Hydrologic model of the saturated zone beneath the Pajarito Plateau using the [DPFEHM](https://github.com/OrchardLANL/DPFEHM.jl) simulator. The model is differentiable and gradients can be efficiently computed with [Zygote](https://github.com/FluxML/Zygote.jl).

To run from in [Julia](https://julialang.org), simply do
```julia
include("ex.jl")
```
The first time it is run it will download a mesh, boundary conditions, etc. from the [waffle2017](https://gitlab.com/LANL-EM/waffle2017) model. Then it runs a steady-state groundwater flow simulation.
