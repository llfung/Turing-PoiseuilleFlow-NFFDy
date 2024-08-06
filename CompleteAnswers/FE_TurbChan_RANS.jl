## Initialise the environment
using Turing, ReverseDiff, StatsPlots, Distributions, Optim, LinearAlgebra, NonlinearSolve
import StatsBase
using MAT

# Plot Package (use PlotlyJS backend)
using Plots
plotlyjs()

# FEM routines
include("../DiscTools/FEM.jl")
using DataInterpolations

# Prandtl Mixing Length Model
    # Based on 
    #    Prandtl, L. (1925). "7. Bericht über Untersuchungen zur ausgebildeten Turbulenz". Z. Angew. Math. Mech. 5 (1): 136–139.
function PrandtlMixingLength(yp,du,p)
    return -(p.κ)^2.0*yp.^2.0 .*abs.(du) .*(du)
end
# Prandtl Mixing Length Model with modification to viscous sublayer
    # Based on 
    #    Prandtl, L. (1925). "7. Bericht über Untersuchungen zur ausgebildeten Turbulenz". Z. Angew. Math. Mech. 5 (1): 136–139.
    #    Van Driest E. (1956), On Turbulent Flow Near a Wall, Journal of Aeronautical Sciences, 1007-1011 
    #    TCACIUC et al. (2022), Bul. Inst. Polit. Iaşi, Vol. 68 (72), Nr. 1 . DOI:10.2478/bipcm-2022-0008 
function PrandtlMixingLengthMod(yp,du,p)
    return -(p.κ)^2.0*yp.^2.0 .*(1.0 .-exp.(-yp/p.Ap)).^2 .*abs.(du) .*(du)
end

function TurbulentChannel_RANS(uavg,p)
    Duavg = FELinDuEst(uavg,dy)
    Duavg[end] = 0.0
    uv = PrandtlMixingLengthMod(ymesh*p.Re,Duavg./p.Re,p)
    # uv = PrandtlMixingLength(ymesh*p.Re,D1*uavg/p.Re,p)
    eqn = D2*uavg/p.Re .- D1*uv .- fmat*ones(Ny)*p.G
    eqn[1] = uavg[1]
    return eqn
end

## Bound of parameters (to prevent blow up and define prior)
logAp_lb= -5.0
logAp_ub= 5.0
logκ_lb = -6.0
logκ_ub = 2.0

## Opening DNS data
DNSData = matread("./Data/LeeMoser2015Channel0180Data.mat");
datamesh = vec(DNSData["y"])
data = vec(DNSData["U"])

## Mesh Set up
Ny = 51 # Global
ymesh,D2,fmat,D1 = FELin(Ny; start = 0.0, stop = 1.0)
dy = 1.0/(Ny-1)

## Parameters
p = (G=-1.0,Re=DNSData["Re"],Ap=exp((logAp_lb+logAp_ub)/2),κ=exp((logκ_lb+logκ_ub)/2))

## Solving the Forward Problem
u0 = (ymesh.-1.0).^2.0 *p.Re/2.0
prob = NonlinearProblem(TurbulentChannel_RANS,u0,p)
sol = solve(prob,SimpleNewtonRaphson();maxiters=100,reltol=1e-8)

## Plotting the solution
plt = plot(DNSData["y"],DNSData["U"],label="DNS",legend=:bottomright)
plot!(plt,ymesh,sol, label="Prandtl Mixing Length (Prior)")

## Data Assimilation
@model function fit(data)
    # LogUniform priors (scale-uniform)
    logAp ~ Uniform(logAp_lb,logAp_ub)
    logκ  ~ Uniform(logκ_lb,logκ_ub)

    # Invoke the solver
    p_local = (G=p.G,Re=p.Re,Ap=exp(logAp),κ=exp(logκ))
    # p_local = (G=p.G,Re=p.Re,κ=exp(logκ))
    prob_local = NonlinearProblem(TurbulentChannel_RANS,u0,p_local)
    sol_local = solve(prob_local,SimpleNewtonRaphson();maxiters=100,reltol=1e-6)
    # Interpolation
    SolInterp_loc = LinearInterpolation(sol_local, ymesh)
    prediction = SolInterp_loc.(datamesh) # piecewise linear interpolation

    # Model how the data is generated
    for i in 2:length(datamesh)
        data[i] ~ Normal(prediction[i], DNSData["U_std"][i])
    end
end

model = fit(data)

## MAP estimation (Only works sometimes)
@time map_estimate= maximum_a_posteriori(model; lb=[logAp_lb,logκ_lb], ub=[logAp_ub, logκ_ub])
StatsBase.coeftable(map_estimate)

## Post-Processing Plots
p_result = (G=p.G,Re=p.Re,Ap=exp(map_estimate.values[:logAp]),κ=exp(map_estimate.values[:logκ]))

prob_result = NonlinearProblem(TurbulentChannel_RANS,u0,p_result)
sol_result = solve(prob_result,SimpleNewtonRaphson();maxiters=100,reltol=1e-8)
plot!(ymesh,sol_result, label="Prandtl Mixing Length (MAP)")