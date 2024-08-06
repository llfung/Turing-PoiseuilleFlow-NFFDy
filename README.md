# Data assimilation workshop with `Turing.jl`
NFFDy workshop on data assimilation with `Turing.jl`, using Poiseuille/Channel flow as example.

## Set up
1. Set up VSCode and the `JULIA` extension.
2. Download the whole repo.
3. Go to the folder and run the following command in `JULIA` before the workshop to set up your environment and get the necessary packages.
```JULIA
] activate VirtualEnv; ] add Turing, StatsPlots, StatsBase, Distributions, Optim, SparseArrays, LinearAlgebra, Plots, PlotlyJS, NonlinearSolve, MAT, DataInterpolations, ForwardDiff, ReverseDiff;  
```
4. Open the folder in VSCode and change the environment to `VirtualEnv` (select the folder `VirtualEnv`).
5. Have fun playing with the scripts!

## Folder Structure
- Data
  - Contain some of the data you'll assimilate, and methods to generate them.
  - You don't have to touch this folder to complete the workshop, but it's fun to play with if you want to explore further.
- DiscTools
  - Contain auxiliary functions, mostly on discretisation. DO NOT TOUCH.
- CompleteAnswers
  - Try not to look into it!
  - The main script is based on Chebyshev Collocation method. There're also Finite Elements version of the same exercise stored in the folder, which you may play with after the workshop.

## Exercises
1. `Cheb_Channel1D_Pressured.jl`
   - A Poiseuille flow flowing through the channel with an unknown pressure gradient `G`. Can we learn about the pressure gradient given some (noisy) meansurements of the flow?
2. `Cheb_Channel1D_Forced.jl`
   - Extending the previous exercise slightly, where instead of a constant pressure gradient across the channel, we now have some unknown arbitrary forcing across the channel. Can we learn about the forcing profile?
3. `Cheb_TurbChan_RANS.jl`
   - Provided some time averaged data of a turbulent channel flow, can we infer the parameters of a modified Prandtl mixing length model for the Reynolds stresses?

## List of packages used
- Turing
- StatsPlots
- StatsBase
- Distributions
- Optim
- SparseArrays
- LinearAlgebra
- Plots
- PlotlyJS
- NonlinearSolve
- MAT
- DataInterpolations
- ForwardDiff
- ReverseDiff