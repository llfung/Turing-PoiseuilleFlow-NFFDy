## Initialise the environment
using Turing, StatsPlots, Distributions, Optim, LinearAlgebra
import StatsBase
using MAT

# Plot Package (use PlotlyJS backend)
using Plots
plotlyjs()

include("./DiscTools/chebdif.jl")

#################### Setting up Forward Problem ####################
## Functions to solve the Poisson Equation
function Poi_Forced(p)
    # LHS: d2/dy^2
    LHS = ## (FILL IN THE BLANK) ##

    # RHS: forcing f
    RHS = ## (FILL IN THE BLANK) ##

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

## Opening data
DataFile = matread("./Data/Channel1D_Forced.mat");

## Interpolation
datamesh = vec(DataFile["y"])
data = vec(DataFile["U"])
ϵ = DataFile["U_std"]

## Plotting the solution
plot(datamesh,data;
    xlabel = "y", 
    ylabel = "U", 
    label="Data",
    legend=:bottomright,
    seriestype=:scatter)

#################### Inverse Problem ####################
## Setting up the model 
@model function fit(data)
    hyperpriorf2 ~ InverseGamma(2,3)

    # Define the prior for the unknown forcing
    f = Vector{Float64}(undef, Ny-2)
    for i in 1:Ny-2
        f[i] ~ Normal(0.0, sqrt(hyperpriorf2))
    end

    # Solve the Foward Problem
    ## (FILL IN THE BLANK) ##
    
    # Interpolation onto datamesh
    ## (FILL IN THE BLANK) ##

    # Model how the data is generated: model solution gives the mean, and noise distribution is Gaussian
    ## (FILL IN THE BLANK) ##
end

model = fit(data)

## MAP estimation
@time map_estimate = maximum_a_posteriori(model)
@show StatsBase.coef(map_estimate)

## Sample chains with forward-mode automatic differentiation (the default).
# @show chain = sample(model, NUTS(), 100; progress=false)
# chain_summary = Turing.summarystats(chain)

## Post-Processing Plots
f_assim = Vector{Float64}(undef, Ny-2)
for i in 1:Ny-2
    f_assim[i] = map_estimate.values[Symbol("f[$i]")]
    # f_assim[i] = chain_summary[Symbol("f[$i]"), :mean]
end

plot(ymesh[2:Ny-1],f_assim; label="Learned f(y)")
plot!(datamesh,vec(DataFile["f"]); label="True f(y)")