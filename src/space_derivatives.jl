using StaticArrays
using LinearAlgebra

@doc raw"""
    posUpwind!(flux, U)

Updates the vector `flux` so that `flux[i]` represents the positive Upwind
flux at U_{i+1/2}, i.e. an Upwind approximation of the flux going into the
left cell from the right cell.

This approximation is only applied for i from 1 to N-1, as the rightmost
flux depends on the boundary condition.
"""
function posUpwind!(flux, U)
    N = length(U)
    @inbounds for i in 1 : N-1
        flux[i] = posUpwind(U, i)
    end
    nothing
end

@doc raw"""
    negUpwind!(flux, U)

Updates the vector `flux` so that `flux[i]` represents the negative Upwind 
flux at U_{i+1/2}, i.e. an Upwind approximation of the flux leaving the
left cell for the right cell.

This approximation is only applied for i from 2 to N, as the leftmost flux
depends on the boundary condition.
"""
function negUpwind!(flux, U)
    N = length(U)
    @inbounds for i in 2 : N
        flux[i] = negUpwind(U, i)
    end
    nothing
end

@doc raw"""
    posWENO3!(flux, U)

Updates the vector "flux" so that flux[i] represents the positive WENO3
flux at U_{i+1/2}, i.e. a WENO3 approximation of the flux going into the
left cell to the right cell.

The rightmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second-to-rightmost flux,
an upwind approximation is used instead.
"""
function posWENO3!(flux, U)
    N = length(U)
    flux[N-1] = posUpwind(U, N-1)
    @inbounds for i in 1 : N-2
        flux[i] = posWENO3(U, i)
    end
    nothing
end

@doc raw"""
    negWENO3!(flux, U)

Updates the vector "flux" so that flux[i] represents the negative WENO3
flux at U_{i+1/2}, i.e. a WENO3 approximation of the flux leaving the left
cell for the right cell.

The rightmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second-to-rightmost flux,
an upwind approximation is used instead.
"""
function negWENO3!(flux, U)
    N = length(U)
    flux[2] = negUpwind(U, 2)
    @inbounds for i in 3 : N
        flux[i] = negWENO3(U, i)
    end
    nothing
end

@doc raw"""
    posWENO5!(flux, U)

Updates the vector "flux" so that flux[i] represents the positive WENO5
flux at U_{i+1/2}, i.e. a WENO5 approximation of the flux going into the
left cell to the right cell.

The rightmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second- and
third-to-rightmost flux, Upwind and WENO3 approximations are used instead.
"""
function posWENO5!(flux, U)
    N = length(U)
    flux[N-1] = posUpwind(U, N-1)
    flux[N-2] = posWENO3(U, N-2)
    @inbounds for i in 1 : N-3
        flux[i] = posWENO5(U, i)
    end
    nothing
end

@doc raw"""
    negWENO5!(flux, U)

Updates the vector "flux" so that flux[i] represents the negative WENO5
flux at U_{i+1/2}, i.e. a WENO5 approximation of the flux leaving the left
cell for the right cell.

The leftmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second- and
third-to-left flux, Upwind and WENO3 approximations are used instead.
"""
function negWENO5!(flux, U)
    N = length(U)
    flux[2] = posUpwind(U, N-1)
    flux[3] = posWENO3(U, N-2)
    @inbounds for i in 4 : N
        flux[i] = negWENO5(U, i)
    end
    nothing
end



@doc raw"""
    posUpwind(U,i::Int)

Computes the positive Upwind flux U_{i+1/2}.
"""
posUpwind(U,i::Int) = U[i+1]

@doc raw"""
    negUpwind(U,i::Int)

Computes the negative Upwind flux U_{i+1/2}.
"""
negUpwind(U,i::Int) = U[i]


@doc raw"""
    posWENO3(U,i::Int)

Computes the positive WENO5 flux U_{i+1/2}.

Ref: https://www.cs.usask.ca/faculty/spiteri/63786-gg.pdf
"""
function posWENO3(U,i::Int)
    local ε = 10^-8

    ## Smoothness indicators
    local β0 = (U[i] - U[i-1])^2
    local β1 = (U[i+1] - U[i])^2

    ## Non-normalized weights
    local α0 = 0.333333333333333333/(ε + β0)^2
    local α1 = 0.666666666666666667/(ε + β1)^2
    local Σα = α0 + α1

    ## Weights (normalized)
    local ω0, ω1 = α0/Σα, α1/Σα

    ## Reconstruction of pointwise values
    local U0 = -0.5*U[i-1] + 1.5*U[i]
    local U1 =  0.5*U[i]   + 0.5*U[i+1]

    ## Final computation
    ω0*U0 + ω1*U1
end

@doc raw"""
    negWENO3(U,i::Int)

Computes the negative WENO5 flux U_{i+1/2}. This could also be obtained by
symmetry from posWENO3.

Ref: https://www.cs.usask.ca/faculty/spiteri/63786-gg.pdf
"""
function negWENO3(U,i::Int)
    local ε = 10^-8

    ## Smoothness indicators
    local β0 = (U[i+2] - U[i+1])^2
    local β1 = (U[i+1] - U[i])^2

    ## Non-normalized weights
    local α0 = 0.333333333333333333/(ε + β0)^2
    local α1 = 0.666666666666666667/(ε + β1)^2
    local Σα = α0 + α1

    ## Weights (normalized)
    local ω0, ω1 = α0/Σα, α1/Σα

    ## Reconstruction of pointwise values
    local U0 = -0.5*U[i+2] + 1.5*U[i+1]
    local U1 =  0.5*U[i+1] + 0.5*U[i]

    ## Final computation
    ω0*U0 + ω1*U1
end

@doc raw"""
    posWENO5(U,i::Int)

Computes the positive WENO5 flux U_{i+1/2}.

Ref: https://www.cs.usask.ca/faculty/spiteri/63786-gg.pdf
"""
function posWENO5(U,i::Int)
    local ε = 10^-8

    ## Smoothness indicators
    local γ = 1.08333333333333333333 # ≈ 13/12
    local β0 = γ*(U[i-2] - 2.0*U[i-1] + U[i])^2 +
                0.25*(U[i-2] - 4.0*U[i-1] + 3.0*U[i])^2
    local β1 = γ*(U[i-1] - 2.0*U[i] + U[i+1])^2 +
                0.25*(U[i-1] - U[i+1])^2
    local β2 = γ*(U[i] - 2.0*U[i+1] + U[i+2])^2 +
                0.25*(3.0*U[i] - 4.0*U[i+1] + U[i+2])^2

    ## Non-normalized weights
    local α0 = 0.3/(ε + β0)^2
    local α1 = 0.6/(ε + β1)^2
    local α2 = 0.1/(ε + β2)^2
    local Σα = α0 + α1 + α2

    ## Weights (normalized)
    local ω0, ω1, ω2 = α0/Σα, α1/Σα, α2/Σα

    ## Reconstruction of pointwise values
    local c0 = SA[ 0.333333333333333333, -1.16666666666666667,  1.833333333333333333]
    local c1 = SA[-0.166666666666666667,  0.83333333333333333,  0.333333333333333333]
    local c2 = SA[ 0.333333333333333333,  0.83333333333333333, -0.166666666666666667]

    local U0 = c0[1]*U[i-2] + c0[2]*U[i-1] + c0[3]*U[i]
    local U1 = c1[1]*U[i-1] + c1[2]*U[i]   + c1[3]*U[i+1]
    local U2 = c2[1]*U[i]   + c2[2]*U[i+1] + c2[3]*U[i+2]

    ## Final computation
    ω0*U0 + ω1*U1 + ω2*U2
end


@doc raw"""
    negWENO5(U,i::Int)

Computes the negative WENO5 flux U_{i+1/2}. This could also be obtained by
symmetry from posWENO5.

Ref: https://www.cs.usask.ca/faculty/spiteri/63786-gg.pdf
"""
function negWENO5(U, i::Int)
    local ε = 10^-8

    ## Smoothness indicators
    local γ = 1.08333333333333333333 # = 13/12
    local β0 = γ*(U[i+1] - 2.0*U[i+2] + U[i+3])^2 +
                0.25*(3.0*U[i+1] - 4.0*U[i+2] + U[i+3])^2
    local β1 = γ*(U[i] - 2.0*U[i+1] + U[i+2])^2 +
                0.25*(U[i] - U[i+2])^2
    local β2 = γ*(U[i-1] - 2.0*U[i] + U[i+1])^2 +
                0.25*(U[i-1] - 4.0*U[i] + 3.0*U[i+1])^2

    ## Non-normalized weights
    local α0 = 0.3/(ε + β0)^2
    local α1 = 0.6/(ε + β1)^2
    local α2 = 0.1/(ε + β2)^2
    local Σα = α0 + α1 + α2

    ## Weights (normalized coefficients)
    local ω0, ω1, ω2 = α0/Σα, α1/Σα, α2/Σα

    ## Reconstructions of pointwise values
    local c0 = SA[ 1.833333333333333333, -1.16666666666666667,  0.333333333333333333]
    local c1  = SA[ 0.333333333333333333,  0.83333333333333333, -0.166666666666666667]
    local c2  = SA[-0.166666666666666667,  0.83333333333333333,  0.333333333333333333]

    local U0 = dot(c0, U[i+1 : i+3])
    local U1 = dot(c1, U[i   : i+2])
    local U2 = dot(c2, U[i-1 : i+1])

    ## Final computation
    ω0*U0 + ω1*U1 + ω2*U2
end