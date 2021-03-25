using LinearAlgebra
using StaticArrays
using OffsetArrays
using Plots
using StiffKinetic

function advection!(dU, flux, U, dx)
    negWENO3!(flux, -U)
    flux[end] = flux[end-1]
    
    @. dU = (flux[1:end] - flux[0:end-1])/dx
end

function RK4(dU, flux, U, dx, dt)
    k1 = -advection!(dU, flux, U, dx)
    k2 = -advection!(dU, flux, U + 0.5*dt*k1, dx)
    k3 = -advection!(dU, flux, U + 0.5*dt*k2, dx)
    k4 = -advection!(dU, flux, U + dt*k3, dx)
    U + dt/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)
end

Nx = 2^7 + 1
x  = range(0.0, 1.0, length=Nx)
dx = step(x)

U = Float64.(0.25 .<= x .<= 0.75)
dU = zero(U)
flux = OffsetVector(zeros(Nx+1), 0:Nx)

T = 1.0
dt = 2^-9
dt_anim = 2^-5
nt_anim = Int(T/dt_anim)

nt_pf = Int(dt_anim/dt)

anim = @animate for t in 0 : dt_anim : T
    plot(x, U, title = "t = $t", ylims=(-0.1, 1.1))
    for _ in 1 : nt_pf
        U .= RK4(dU, flux, U, dx, dt)
    end
end

gif(anim, "advection.gif", fps=10)