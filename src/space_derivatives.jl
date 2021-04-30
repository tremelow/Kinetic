using StaticArrays
using OffsetArrays
import Base.Iterators: drop, take
import LinearAlgebra: dot

include("space_derivatives/upwind.jl")
include("space_derivatives/weno3.jl")
include("space_derivatives/weno5.jl")
include("space_derivatives/second_order.jl")


function compute_flux!(U :: ScalarQuantity)
    posWENO3!(U.posFlux, U.data)
    negWENO3!(U.negFlux, U.data)
    U.posFlux[0], U.negFlux[end] = U.posFlux[1], U.negFlux[end-1]
end

function compute_fluxP!(pb :: RelaxPb)
    compute_flux!(pb.U)
    @. pb.fluxP = 0.5 * (pb.p(pb.U.posFlux) + pb.p(pb.U.negFlux))
    nothing
end

function compute_fluxV!(pb :: RelaxPb)
    compute_flux!(pb.V)
    @. pb.fluxV = 0.5 * (pb.V.posFlux + pb.V.negFlux)
    nothing
end


function compute_dxP!(pb :: RelaxPb, dx⁻¹)
    compute_flux!(pb.U)
    @. pb.fluxP = 0.5 * ( pb.p(pb.U.posFlux) + pb.p(pb.U.negFlux) )
    
    for i in 1 : length(pb.dxP)
        pb.dxP[i] = (pb.fluxP[i] - pb.fluxP[i-1]) * dx⁻¹
    end
    
    nothing
end

function compute_dxV!(pb :: RelaxPb, dx⁻¹)
    compute_flux!(pb.V)
    @. pb.fluxV = 0.5 * ( pb.V.posFlux + pb.V.negFlux )
    
    for i in 1 : length(pb.dxV)
        pb.dxV[i] = (pb.fluxV[i] - pb.fluxV[i-1]) * dx⁻¹
    end

    nothing
end
