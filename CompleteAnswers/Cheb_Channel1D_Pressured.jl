## Initialise the environment
using Turing, StatsPlots, Distributions, Optim
import StatsBase

# Plot Package (use PlotlyJS backend)
using Plots
plotlyjs()

include("../DiscTools/chebdif.jl")

## Functions to solve the Poisson Equation
function Poi(p)
    # LHS: d2/dy^2
    LHS = [ view(D1,1,1:Ny)'; # BC at y=1 (Neumann)
            view(D2,2:Ny-1,1:Ny);
            zeros(Ny-1)' [1.0] ] # BC at y=0 (no slip)

    # RHS: Re*dp/dx (=G)
    RHS = [ 0.0; # BC at y=1 (Neumann)
            p.G*ones(Ny-2);
            0.0 ] # BC at y=0 (no slip)

    sol = LHS\RHS # Direct inversion, linear problem

    return sol
end

## Setting up and solving the Forward Problem
# Initialise the Mesh - y ∈ [0,1] (half channel)
Ny = 21 # Global
ymesh,DM = chebdif(Ny,2)
ymesh = (ymesh.+1.0)./2.0
D1 = DM[:,:,1]*2.0 # Global
D2 = DM[:,:,2]*4.0 # Global

# Parameters
p = (G=-10.0,)

# Solving the Poisson problem
sol = Poi(p)

# Randomly sample data points from the forward solution
NData = 10
datamesh = rand(NData) # Random sampling location in y ∈ [0,1]
sol_int = chebint(sol,(datamesh.*2.0.-1.0)) # Interpolate the forward solution from the random sampling location (datamesh need to be mapped from [0,1] to [-1,1] for the chebyshev method to work)

# Add noise to the interpolated solution to get the data
ϵ = 0.1
data = sol_int .+ ϵ*randn(NData) 

## Plotting the solution
plot(ymesh,sol;
    xlabel = "y", 
    ylabel = "U", 
    label="Model",
    linewidth=2,
    color=:black,
    legend=:bottomright,)
plot!(datamesh,data;
    label="Data", 
    seriestype=:scatter)

## Data Assimilation
@model function fit(data)
    # OPTION 1: Define variance as hyperpriors
    # σ2 ~ InverseGamma(2, 3)
    # GScale ~ InverseGamma(2,3)

    # OPTION 2: Define variance by known noise and scaling argument 
    # noise level given by ϵ
    σ2 = ϵ^2.0
    # Scale of G given by U / L^2
    GScale2 = (maximum(abs,data) / (1.0^2.0))^2.0

    # Define the prior for the unknown pressure gradient
    G ~ Normal(0.0, sqrt(GScale2))

    # Solve the Poisson problem with the parameter G
    p = (G=G,)
    sol = Poi(p)
    
    # Interpolation onto datamesh
    sol_int = chebint(sol,(datamesh.*2.0.-1.0))

    # Model how the data is generated: model solution gives the mean, and noise distribution is Gaussian
    data ~ MvNormal(sol_int, σ2*I)
    # Alternative way to write the above line
    # for i in 1:NData
    #     data[i] ~ Normal(sol_int[i], sqrt(σ2))
    # end
end

model = fit(data)

## MLE estimation
# @time mle_estimate = maximum_likelihood(model)
# StatsBase.coeftable(mle_estimate)

## MAP estimation
@time map_estimate = maximum_a_posteriori(model)
StatsBase.coeftable(map_estimate)

## Sample chains with forward-mode automatic differentiation (the default).
# @show chain = sample(model, NUTS(), 1000; progress=false)
# plot(chain)