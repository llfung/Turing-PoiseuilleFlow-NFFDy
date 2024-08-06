## Stencil for piecewise linear finite element method in uniform 1D mesh 
using SparseArrays

function FELin(Nx; start = -1.0, stop = 1.0)
    xmesh = collect(range(start, stop; length=Nx))
    dx = (stop-start)/(Nx-1)
    lap_mat = spdiagm(Nx,Nx,-1=>ones(Nx-1), 0=>-2*[0.5;ones(Nx-2);0.5], 1=>ones(Nx-1))/dx^2
    forcing_mat = spdiagm(Nx,Nx,-1=>(1.0/6.0)*ones(Nx-1), 0=>(2.0/3.0)*[0.5;ones(Nx-2);0.5], 1=>(1.0/6.0)*ones(Nx-1))
    D1_mat = spdiagm(Nx,Nx,-1=>-ones(Nx-1), 0=>[-1.0;zeros(Nx-2);1.0], 1=>ones(Nx-1))/2.0/dx
    return xmesh,lap_mat,forcing_mat,D1_mat
end

function FELinDuEst(u,dx)
    du=(u[3:end]-u[1:end-2])/2/dx
    return [(u[2]-u[1])/dx;du;(u[end]-u[end-1])/dx]
end