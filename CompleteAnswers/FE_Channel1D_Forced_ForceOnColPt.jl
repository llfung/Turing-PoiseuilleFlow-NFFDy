## Initialise the environment
using Turing, StatsPlots, Distributions, Optim, LinearAlgebra
import StatsBase

# Plot Package (use PlotlyJS backend)
using Plots
plotlyjs()

include("../DiscTools/FEM.jl")
using DataInterpolations

## Functions to solve the Poisson Equation
function Poi_Forced(p)
    # LHS: d2/dy^2
    LHS = view(lap,2:Ny,2:Ny);

    # RHS: forcing f
    RHS_Full = fmat*p.f;
    RHS = view(RHS_Full,2:Ny)

    sol = LHS\RHS # Direct inversion, linear problem

    return [0.0;sol]
end

## Setting up and solving the Forward Problem
# Initialise the Mesh - y ∈ [0,1] (half channel)
Ny = 21 # Global
ymesh,lap,fmat = FELin(Ny; start = -0.0, stop = 1.0)

# Parameters (forcing function)
f = (4*(ymesh .-1.0).^2.0 .-3*(ymesh .-1.0).^4.0 .-0.5)*10.0
# Solving the PDE
sol = Poi_Forced((f=f,))

# Add noise to the solution
ϵ = 0.01
data = zeros(Ny)
data[2:Ny] = sol[2:Ny] .+ ϵ*randn(Ny-1)

## Plotting the solution
plot(ymesh,sol;
    xlabel = "y", 
    ylabel = "U", 
    label="Model",
    linewidth=2,
    color=:black,
    legend=:bottomright,)
plot!(ymesh,data;
    label="Data", 
    seriestype=:scatter)

## Data Assimilation
@model function fit(data)
    # σ2 ~ InverseGamma(2, 3)
    σ2 = ϵ^2.0

    hyperpriorf2 ~ InverseGamma(2,3)

    # Define the prior for the unknown forcing
    f = Vector{Float64}(undef, Ny)
    for i in 1:Ny
        f[i] ~ Normal(0.0, sqrt(hyperpriorf2))
    end

    # Solve the Poisson problem with the parameter G
    p = (f=f,)
    sol = Poi_Forced(p)

    # Interpolation onto datamesh
    # SolInterp = LinearInterpolation(sol, ymesh)
    # sol_int = SolInterp.(datamesh) # piecewise linear interpolation

    # Model how the data is generated: model solution gives the mean, and noise distribution is Gaussian
    # data ~ MvNormal(sol_int, σ2*I)
    # Alternative way to write the above line
    for i in 2:length(data)
        data[i] ~ Normal(sol[i], sqrt(σ2))
    end
end

model = fit(data)

## MAP estimation
@time map_estimate = maximum_a_posteriori(model)
@show map_estimate.lp
@show StatsBase.coef(map_estimate)
# StatsBase.vcov(map_estimate)
# StatsBase.informationmatrix(map_estimate)
# StatsBase.coeftable(map_estimate)

## Sample chains with forward-mode automatic differentiation (the default).
# @show chain = sample(model, NUTS(), 400; progress=false)
# chain_summary = Turing.summarystats(chain)
# plot(chain)

## Post-Processing Plots - Extracting and plotting the mean of the forcing function
f_assim = Vector{Float64}(undef, Ny)
for i in 1:Ny
    f_assim[i] = map_estimate.values[Symbol("f[$i]")]
    # f_assim[i] = chain_summary[Symbol("f[$i]"), :mean]
end

plot(ymesh,f_assim)
plot!(ymesh,f)