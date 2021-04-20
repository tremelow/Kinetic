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

    dxU :: Vector
    dxV :: Vector

    posFluxDxU :: OffsetVector
    negFluxDxU :: OffsetVector
    rosFluxDxU :: OffsetVector
    posFluxDxV :: OffsetVector
    negFluxDxV :: OffsetVector
    rosFluxDxV :: OffsetVector

    dx2_U :: Vector

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

        posFluxDxU = OffsetVector(zeros(Nx+1), 0:Nx)
        negFluxDxU = OffsetVector(zeros(Nx+1), 0:Nx)
        rosFluxDxU = OffsetVector(zeros(Nx+1), 0:Nx)
        posFluxDxV = OffsetVector(zeros(Nx+1), 0:Nx)
        negFluxDxV = OffsetVector(zeros(Nx+1), 0:Nx)
        rosFluxDxV = OffsetVector(zeros(Nx+1), 0:Nx)

        dU, dV = zero(U), zero(V)
        dxU, dxV = zero(U), zero(V)
        dx2_U = zero(U)


        new(copy(U),copy(V),
            posFluxU,negFluxU,rosFluxU,
            posFluxV,negFluxV,rosFluxU,
            dxU,dxV,
            posFluxDxU,negFluxDxU,rosFluxDxU,
            posFluxDxV,negFluxDxV,rosFluxDxU,
            dx2_U,dU,dV)
    end
end

Nx = 2^10 + 1
x  = range(0.0, 4.0, length=Nx)
dx = step(x)

uL, uR = 4.0, 2.0
vL, vR = 0.0, 0.0

U = uL*Float64.(0.0 .<= x .<= 2.0) + uR*Float64.(2.0 .< x .<= 4.0)
V = vL*Float64.(0.0 .<= x .<= 2.0) + vR*Float64.(2.0 .< x .<= 4.0)
pb = RelaxPb(U,V)

ε = 2^-20

function imex_bdf1_du!(pb :: RelaxPb, dt, dx, ε)
    posUpwind!(pb.posFluxU, pb.U)
    posUpwind!(pb.posFluxV, pb.V)
    negUpwind!(pb.negFluxU, pb.U)
    negUpwind!(pb.negFluxV, pb.V)
    pb.posFluxU[0], pb.negFluxU[end] = uL, uR
    pb.posFluxV[0], pb.negFluxV[end] = vL, vR

    
    local Θ = 1/(ε^2 + dt)
    local rosΘ = 0.0
    @. pb.rosFluxU = pb.posFluxU + pb.negFluxU +
                    -rosΘ * (pb.posFluxV - pb.negFluxV)
    @. pb.rosFluxV = pb.posFluxV + pb.negFluxV +
                    -rosΘ * (pb.posFluxU - pb.negFluxU)

    local dx⁻¹ = 1.0 / dx
    pb.dxU .= 0.5 * dx⁻¹ * diff(pb.rosFluxU.parent)
    pb.dxV .= 0.5 * dx⁻¹ * diff(pb.rosFluxV.parent)


    posUpwind!(pb.posFluxDxU, pb.dxU)
    posUpwind!(pb.posFluxDxV, pb.dxV)
    negUpwind!(pb.negFluxDxU, pb.dxU)
    negUpwind!(pb.negFluxDxV, pb.dxV)
    pb.posFluxDxU[0], pb.negFluxDxU[end] = 0.0, 0.0
    pb.posFluxDxV[0], pb.negFluxDxV[end] = 0.0, 0.0

    @. pb.rosFluxDxU = pb.posFluxDxU + pb.negFluxDxU +
                    -rosΘ * (pb.posFluxDxV - pb.negFluxDxV)

    pb.dx2_U .= 0.5 * dx⁻¹ * diff(pb.rosFluxDxU.parent)

    # fd_diff2_ord6!(pb.dx2_U, U, dx⁻¹)

    @. pb.dU =  dt * Θ * (ε^2 * pb.dxV + dt * pb.dx2_U)
    @. pb.dV = -dt * Θ * ( pb.V + pb.dxU )

    pb.dU[1], pb.dU[end] = 0.0, 0.0
    pb.dV[1], pb.dV[end] = 0.0, 0.0

    nothing
end

T = 0.25
dt = (2^-5) * dx^2

imex_bdf1_du!(pb, dt, dx, ε)

dt_anim = 2^-6 * T
nt_anim = Int(T/dt_anim)

nt_pf = Int(dt_anim/dt)

anim = @animate for t in 0 : dt_anim : T
    plot(x, pb.U, title = "t = $t", ylims=(1.5, 5.5))
    for _ in 1 : nt_pf
        imex_bdf1_du!(pb, dt, dx, ε)
        @. pb.U += pb.dU
        @. pb.V += pb.dV
    end
end

gif(anim, "riemann.gif", fps=10)
