###############################
#
#           UPWIND
#
###############################


##-----------------------
## Positive

@doc raw"""
    posUpwind!(flux, U)

Updates the vector `flux` so that `flux[i]` represents the positive Upwind
flux at x_{i+1/2}, i.e. an Upwind approximation of the flux going into the
right cell from the left cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

This approximation is only applied for i from 1 to N, as the leftmost
flux (at x_{1/2}) depends on the boundary condition.
"""
function posUpwind!(flux, U)
    @. flux[1 : end] = U
    nothing
end

@doc raw"""
    posUpwind(U,i::Int)

Computes the positive Upwind flux x_{i+1/2}.
"""
posUpwind(U,i::Int) = U[i]


##-----------------------
## Negative

@doc raw"""
    negUpwind!(flux, U)

Updates the vector `flux` so that `flux[i]` represents the negative Upwind
flux at x_{i+1/2}, i.e. an Upwind approximation of the flux going from the
right cell to the left cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

This approximation is only applied for i from 1 to N-1, as the rightmost
flux depends on the boundary condition.
"""
function negUpwind!(flux, U)
    @. flux[0 : end-1] = U
    nothing
end

@doc raw"""
    negUpwind(U,i::Int)

Computes the negative Upwind flux x_{i+1/2}.
"""
negUpwind(U, i::Int) = U[i+1]