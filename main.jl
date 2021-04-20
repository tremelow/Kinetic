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

struct ScalarQuantity
    data :: Vector
    posFlux :: OffsetVector
    negFlux :: OffsetVector

    function ScalarQuantity(U)
        Nx = length(U)
        posFlux = OffsetVector(zeros(Nx+1), 0:Nx)
        negFlux = OffsetVector(zeros(Nx+1), 0:Nx)

        new(copy(U), posFlux, negFlux)
    end
end

function compute_flux!(U :: ScalarQuantity, uL, uR)
    posWENO5!(U.posFlux, U.data)
    negWENO5!(U.negFlux, U.data)
    U.posFlux[0], U.negFlux[end] = uL, uR
end

struct DoubleQuantity
    U :: ScalarQuantity
    V :: ScalarQuantity
    rosFluxU :: OffsetVector
    rosFluxV :: OffsetVector

    function DoubleQuantity(U,V)
        Nx = length(U)
        @assert Nx == length(V)
        rosFluxU = OffsetVector(zeros(Nx+1), 0:Nx)
        rosFluxV = OffsetVector(zeros(Nx+1), 0:Nx)

        new(ScalarQuantity(U), ScalarQuantity(V), rosFluxU, rosFluxV)
    end
end

function compute_ros_flux!(uv :: DoubleQuantity, uL, uR, vL, vR, Θ)
    compute_flux!(uv.U, uL, uR)
    compute_flux!(uv.V, vL, vR)
    @. uv.rosFluxU = 0.5 * (uv.U.posFlux + uv.U.negFlux +
                            -Θ * (uv.V.posFlux - uv.V.negFlux) )
    @. uv.rosFluxV = 0.5 * (uv.V.posFlux + uv.V.negFlux +
                            -Θ * (uv.U.posFlux - uv.U.negFlux) )
end


function assign_dx!(dx_uv, uv, uL, uR, vL, vR, rosΘ, dx⁻¹)
    compute_ros_flux!(uv, uL, uR, vL, vR, rosΘ)
    dx_uv.U.data .= dx⁻¹ * diff(pb.uv.rosFluxU.parent)
    dx_uv.V.data .= dx⁻¹ * diff(pb.uv.rosFluxV.parent)
end



struct RelaxPb
    uv :: DoubleQuantity
    dx_uv :: DoubleQuantity
    dx2_uv :: DoubleQuantity

    dU :: Vector
    dV :: Vector

    function RelaxPb(U,V)
        dx_U, dx_V = zero(U), zero(V)
        dx2_U, dx2_V = zero(U), zero(V)
        dU, dV = zero(U), zero(V)

        new(DoubleQuantity(U, V), DoubleQuantity(dx_U, dx_V),
            DoubleQuantity(dx2_U, dx2_V), dU, dV)
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

ε = 2^-1

function imex_bdf1_du!(pb :: RelaxPb, dt, dx, ε)
    local Θ = 1.0 / (ε^2 + dt)
    local dx⁻¹ = 1.0 / dx

    local rosΘ = 0.0
    assign_dx!(pb.dx_uv, pb.uv, uL, uR, vL, vR, rosΘ, dx⁻¹)
    assign_dx!(pb.dx2_uv, pb.dx_uv, 0, 0, 0, 0, rosΘ, dx⁻¹)

    # fd_diff2_ord6!(pb.dx2_U, U, dx⁻¹)

    @. pb.dU =  dt * Θ * (ε^2 * pb.dx_uv.V.data + dt * pb.dx2_uv.U.data)
    @. pb.dV = -dt * Θ * ( pb.uv.V.data + pb.dx_uv.U.data )

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
    plot(x, pb.uv.U.data, title = "t = $t", ylims=(1.5, 5.5))
    for _ in 1 : nt_pf
        imex_bdf1_du!(pb, dt, dx, ε)
        @. pb.uv.U.data += pb.dU
        @. pb.uv.V.data += pb.dV
    end
end

gif(anim, "riemann.gif", fps=10)
