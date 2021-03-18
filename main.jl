using LinearAlgebra
using StaticArrays
using CircularArrays
using Plots
using StiffKinetic

function advection(U, dx)
    flux = copy(U)
    for i in 1 : length(U)
        flux[i] = (posWENO5(U,i-1) - posWENO5(U,i))/dx
    end
    flux
end

function RK4(U, dx, dt)
    k1 = advection(U, dx)
    k2 = advection(CircularVector(U + 0.5*dt*k1), dx)
    k3 = advection(CircularVector(U + 0.5*dt*k2), dx)
    k4 = advection(CircularVector(U + dt*k3), dx)
    CircularVector(U + dt/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4))
end

Nx = 2^7 + 1
x  = range(0.0, 1.0, length=Nx)
dx = step(x)

U = CircularVector(Float64.(0.25 .<= x .<= 0.75))

T = 2.0
dt = 2^-9
dt_anim = 2^-5
nt_anim = Int(T/dt_anim)

nt_pf = Int(dt_anim/dt)

anim = @animate for t in 0 : dt_anim : T
    plot(x, Vector(U), title = "t = $t")
    for _ in 1 : nt_pf
        U .= RK4(U, dx, dt)
    end
end

gif(anim, "advection.gif", fps=10)