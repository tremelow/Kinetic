using LinearAlgebra
using StaticArrays
using OffsetArrays
using Plots
using StiffKinetic
using Printf


Nx = 2^8 + 1
x  = range(0.0, 4.0, length=Nx)
dx = step(x)

uL, uR = 4.0, 2.0
U0 = uL*Float64.(0.0 .<= x .<= 2.0) + uR*Float64.(2.0 .< x .<= 4.0)
V0 = zero(U0)

p = u -> u
f = u -> u # /!\ we may use |f'| = 1 in the schemes

pb = RelaxPb(U0,V0,f,p)
pbs = [pb]


ε = 2^-1

T = 0.25
dt = (2^-4) * dx^2

dt_anim = 2^-6 * T
nt_anim = Int(T/dt_anim)

nt_pf = Int(dt_anim/dt)

anim = @animate for t in 0 : dt_anim : T
    str_t = @sprintf "t = %1.3f" t
    plot(x, pbs[1].U.data, title = str_t, ylims=(1.8, 4.2))
    # plot!(x, pbs[1].V.data, title = "t = $t", ylims=(-5, 4.2))
    for _ in 1 : nt_pf
        imex_bdf4!(pbs, dt, dx, ε)
    end
end

gif(anim, "riemann.gif", fps=24)
