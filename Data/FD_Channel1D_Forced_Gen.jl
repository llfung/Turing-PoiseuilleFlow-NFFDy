using SparseArrays
# Defining differentiation matrix
function D2(Ny,dx)
    return spdiagm(Ny,Ny,-1=>ones(Ny-1), 0=>-2*ones(Ny), 1=>ones(Ny-1))/dx^2
end

# Function to solve 
function Poi(Ny,f; start=-1.0, stop=1.0)
    dx = (stop-start)/(Ny-1)
    lap = D2(Ny,dx)
    LHS = view(lap, 2:Ny-1, 2:Ny-1)
    RHS = view(f,2:Ny-1)

    sol = LHS\RHS
    return [0.0;sol;0.0]
end

## Running the solver
Ny = 51 # ODD
ymesh = collect(range(-1.0, 1.0; length=Ny))

# Parameters (forcing function)
f = (4*ymesh.^2.0 .-3*ymesh.^4.0 .-0.5)*10.0
# Solving the PDE
sol = Poi(Ny,f)

# Add noise to the solution
ϵ = 0.01
data = zeros(Ny)
data[2:Ny-1] = sol[2:Ny-1] .+ ϵ*randn(Ny-2)


## Plotting the solution
using Plots

plot(ymesh,sol;
    xlabel = "y", 
    ylabel = "u", 
    label="Model",
    linewidth=2,
    color=:black,
    legend=:bottomright,)
plot!(ymesh,data;
    label="Data", 
    seriestype=:scatter)

# Save the data
using MAT
datamesh_out = ymesh[2:Int((Ny+1)//2)] .+1.0
data_out = data[2:Int((Ny+1)//2)]
f_out = f[2:Int((Ny+1)//2)]
matopen("Channel1D_Forced.mat", "w") do fil
    write(fil, "y", datamesh_out)
    write(fil, "U", data_out)
    write(fil, "U_std",ϵ)
    write(fil, "f", f_out)
end