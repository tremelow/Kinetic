using LinearAlgebra
using StaticArrays
using OffsetArrays
using Plots
using StiffKinetic


"""
Consider the system
    ∂ₜu + ∂ₓv = 0
    ∂ₜv = -(v + ∂ₓu)/ε²
"""

struct RelaxPb
    U :: Vector
    V :: Vector

    posFluxU :: OffsetVector
    negFluxU :: OffsetVector
    rosFluxU :: OffsetVector
    posFluxV :: OffsetVector
    negFluxV :: OffsetVector
    rosFluxV :: OffsetVector

    dx2_U :: OffsetVector

    dU :: Vector
    dV :: Vector

    function RelaxPb(U,V)
        Nx = length(U)
        @assert Nx == length(V)

        posFluxU = OffsetVector(zeros(Nx+1), 0:Nx)
        negFluxU = OffsetVector(zeros(Nx+1), 0:Nx)
        rosFluxU = OffsetVector(zeros(Nx+1), 0:Nx)
        posFluxV = OffsetVector(zeros(Nx+1), 0:Nx)
        negFluxV = OffsetVector(zeros(Nx+1), 0:Nx)
        rosFluxV = OffsetVector(zeros(Nx+1), 0:Nx)

        dx2_U = zero(U)

        dU, dV = zero(U), zero(V)

        new(copy(U),copy(V),posFluxU,negFluxU,rosFluxU,
            posFluxV,negFluxV,rosFluxU,dx2_U,dU,dV)
    end
end

Nx = 2^7 + 1
x  = range(0.0, 1.0, length=Nx)
dx = step(x)

U = Float64.(0.25 .<= x .<= 0.75)
V = Float64.(0.25 .<= x .<= 0.75)
pb = RelaxPb(U,V)

function imex_bdf1_du(pb :: RelaxPb, dt, dx, ε)
    posWENO5!(pb.posFluxU, U)
    posWENO5!(pb.posFluxV, V)
    negWENO5!(pb.negFluxU, U)
    negWENO5!(pb.negFluxV, V)



    local Θ = 1/(ε^2 + dt)
    local dx⁻¹ = 1.0 / dx
    local ∂ₓu = 0.5 * dx⁻¹ * (pb.posFluxU + pb.negFluxU 
                                - Θ * (pb.posFluxV - pb.negFluxV))
    local ∂ₓv = 0.5 * dx⁻¹ * (pb.posFluxV + pb.negFluxV 
                                - Θ * (pb.posFluxU - pb.negFluxU))
    dV .= -dt/(ε^2 + dt) * ( V + ∂ₓu )

    dU = dt / (ε^2 + dt) * ∂ₓV # with new V
end

T = 1.0
dt = 2^-9
dt_anim = 2^-6
nt_anim = Int(T/dt_anim)

nt_pf = Int(dt_anim/dt)

anim = @animate for t in 0 : dt_anim : T
    plot(x, U, title = "t = $t", ylims=(-0.1, 1.1))
    for _ in 1 : nt_pf
        U .= imex_bdf1_du(dU, dV, fluxU, fluxV, U, V, dt, dx, ε)
    end
end

gif(anim, "advection.gif", fps=10)