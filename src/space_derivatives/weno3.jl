#################################
#
#           WENO3
#
##################################


##-------------------------
## Positive

@doc raw"""
    posWENO3!(flux, U)

Updates the vector "flux" so that flux[i] represents the positive WENO3
flux at x_{i+1/2}, i.e. a WENO3 approximation of the flux going from the
left cell to the right cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

The leftmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second-to-leftmost and the
rightmost fluxes, an Upwind approximation is used instead.
"""
function posWENO3!(flux, U)
    N = length(U)
    flux[1] = posUpwind(U, 1)
    flux[N] = posUpwind(U, N)

    local ε = 10^-8
    # Squared smoothness indicators ("non-zero"-ified)
    local β = diff(U).^2
    @. β = 1.0 / ( ε + β )^2
    # Non-normalized weights
    local ω0 = 0.333333333333333333 .* β[1:N-2]
    local ω1 = 0.666666666666666667 .* β[2:N-1]
    local Σα = ω0 + ω1
    # Normalize weights
    @. ω0 /= Σα
    @. ω1 /= Σα
    # Pointwise values
    local U0 = -0.5*U[1:N-2] + 1.5*U[2:N-1]
    local U1 =  0.5*U[2:N-1] + 0.5*U[3:N]

    @. flux[2:N-1] = ω0*U0 + ω1*U1

    nothing
end


@doc raw"""
    posWENO3(U,i::Int)

Computes the positive WENO5 flux of U at x_{i+1/2}.

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

##---------------------------
## Negative

@doc raw"""
    negWENO3!(flux, U)

Updates the vector "flux" so that flux[i] represents the negative WENO3
flux at x_{i+1/2}, i.e. a WENO3 approximation of the flux going from the
right cell to the left cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

The rightmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second-to-rightmost and the
leftmost fluxes, an upwind approximation is used instead.
"""
function negWENO3!(flux, U)
    N = length(U)
    flux[0] = negUpwind(U, 0)
    flux[N-1] = negUpwind(U, N-1)

    local ε = 10^-8
    # Inverse of squared smoothness indicators ("non-zero"-ified)
    local β = diff(U) .^2
    @. β = 1.0 / (ε + β)^2
    # Non-normalized weights
    local ω0 = 0.333333333333333333 .* β[2:N-1]
    local ω1 = 0.666666666666666667 .* β[1:N-2]
    local Σα = ω0 + ω1
    # Normalize weights
    @. ω0 /= Σα
    @. ω1 /= Σα
    # Pointwise values
    local U0 = -0.5*U[3:N]   + 1.5*U[2:N-1]
    local U1 =  0.5*U[2:N-1] + 0.5*U[1:N-2]

    @. flux[1:N-2] = ω0*U0 + ω1*U1

    # @inbounds for i in 2 : N-2
    #     flux[i] = negWENO3(U, i)
    # end
    nothing
end

@doc raw"""
    negWENO3(U,i::Int)

Computes the negative WENO5 flux of U at x_{i+1/2}. This could also be
obtained by symmetry from posWENO3.

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
