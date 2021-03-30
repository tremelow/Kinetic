###############################
#
#           WENO5
#
###############################


##---------------------------
## Positive

@doc raw"""
    posWENO5!(flux, U)

Updates the vector "flux" so that flux[i] represents the positive WENO5
flux at x_{i+1/2}, i.e. a WENO5 approximation of the flux going from the
left cell to the right cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

The leftmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second- and
third-to-leftmost fluxes as well as the rightmost and second-to-rightmost
fluxes, Upwind and WENO3 approximations are used instead.
"""
function posWENO5!(flux, U)
    N = length(U)
    flux[1] = posUpwind(U, 1)
    flux[2] = posWENO3(U, 2)
    flux[N-1] = posWENO3(U, N-1)
    flux[N]   = posUpwind(U, N)


    local ε = 10^-8
    ## Smoothness indicators
    local S0 = take(zip(U, drop(U,1), drop(U,2)), N-4) # Stencils for β0
    local S1 = take(zip(drop(U,1), drop(U,2), drop(U,3)), N-4)
    local S2 = zip(drop(U,2), drop(U,3), drop(U,4))

    local γ = 1.08333333333333333333 # ≈ 13/12
    
    local β0 = [γ*(Ug - 2.0*Um + Ud)^2 + 0.25*(Ug - 4.0*Um + 3.0*Ud)^2
                for (Ug, Um, Ud) ∈ S0]

    local β1 = [γ*(Ug - 2.0*Um + Ud)^2 + 0.25*(Ug - Ud)^2
                for (Ug, Um, Ud) ∈ S1]

    local β2 = [γ*(Ug - 2.0*Um + Ud)^2 + 0.25*(3.0*Ug - 4.0*Um + Ud)^2
                for (Ug, Um, Ud) ∈ S2]

    ## Non-normalized weights
    local ω0 = 0.3 ./ (ε .+ β0).^2
    local ω1 = 0.6 ./ (ε .+ β1).^2
    local ω2 = 0.1 ./ (ε .+ β2).^2
    # normalize them
    local Σα = ω0 + ω1 + ω2
    @. ω0 /= Σα
    @. ω1 /= Σα
    @. ω2 /= Σα

    ## Pointwise values
    local c0 = SA[ 0.333333333333333333, -1.16666666666666667,  1.833333333333333333]
    local c1 = SA[-0.166666666666666667,  0.83333333333333333,  0.333333333333333333]
    local c2 = SA[ 0.333333333333333333,  0.83333333333333333, -0.166666666666666667]

    local U0 = [c0[1]*Ug + c0[2]*Um + c0[2]*Ud for (Ug, Um, Ud) in S0]
    local U1 = [c1[1]*Ug + c1[2]*Um + c1[2]*Ud for (Ug, Um, Ud) in S1]
    local U2 = [c2[1]*Ug + c2[2]*Um + c2[2]*Ud for (Ug, Um, Ud) in S2]

    @. flux[3:N-2] = ω0*U0 + ω1*U1 + ω2*U2

    nothing
end



@doc raw"""
    posWENO5(U,i::Int)

Computes the positive WENO5 flux of U at x_{i+1/2}.

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


##-----------------------
## Negative

@doc raw"""
    negWENO5!(flux, U)

Updates the vector "flux" so that flux[i] represents the negative WENO5
flux at x_{i+1/2}, i.e. a WENO5 approximation of the flux going from the
left cell to the right cell. Necessarily, `flux` is supposed to be an
OffsetArray with indices from 0 to N.

The rightmost flux is not computed, as it depends on the boundary
condition. For the stencil to be valid for the second- and
third-to-rightmost fluxes as well as the leftmost and second-to-leftmost
fluxes, Upwind and WENO3 approximations are used instead.
"""
function negWENO5!(flux, U)
    N = length(U)
    flux[0] = negUpwind(U, 0)
    flux[1] = negWENO3(U, 1)
    flux[N-2] = negWENO3(U, N-2)
    flux[N-1] = negUpwind(U, N-1)
    
    @inbounds for i in 2 : N-3
        flux[i] = negWENO5(U, i)
    end
    
    nothing
end

@doc raw"""
    negWENO5(U,i::Int)

Computes the negative WENO5 flux of U at x_{i+1/2}. This could also be
obtained by symmetry from posWENO5.

Ref: https://www.cs.usask.ca/faculty/spiteri/63786-gg.pdf
"""
function negWENO5(U, i::Int)
    local ε = 10^-8

    ## Smoothness indicators
    local γ = 1.08333333333333333333 # = 13/12
    local β0 = γ*(U[i+3] - 2.0*U[i+2] + U[i+1])^2 +
                0.25*(U[i+3] - 4.0*U[i+2] + 3.0*U[i+1])^2
    local β1 = γ*(U[i+2] - 2.0*U[i+1] + U[i])^2 +
                0.25*(U[i+2] - U[i])^2
    local β2 = γ*(U[i+1] - 2.0*U[i] + U[i-1])^2 +
                0.25*(3.0*U[i+1] - 4.0*U[i] + U[i-1])^2

    ## Non-normalized weights
    local α0 = 0.3/(ε + β0)^2
    local α1 = 0.6/(ε + β1)^2
    local α2 = 0.1/(ε + β2)^2
    local Σα = α0 + α1 + α2

    ## Weights (normalized coefficients)
    local ω0, ω1, ω2 = α0/Σα, α1/Σα, α2/Σα

    ## Reconstructions of pointwise values
    # local c0 = SA[ 1.833333333333333333, -1.16666666666666667,  0.333333333333333333]
    # local c1  = SA[ 0.333333333333333333,  0.83333333333333333, -0.166666666666666667]
    # local c2  = SA[-0.166666666666666667,  0.83333333333333333,  0.333333333333333333]
    
    # local U0 = c0[1]*U[i+1] + c0[2]*U[i+2] + c0[3]*U[i+3]
    # local U1 = c1[1]*U[i]   + c1[2]*U[i+1] + c1[3]*U[i+2]
    # local U2 = c2[1]*U[i-1] + c2[2]*U[i]   + c2[3]*U[i+1]

    local c0 = SA[ 0.333333333333333333, -1.16666666666666667,  1.833333333333333333]
    local c1 = SA[-0.166666666666666667,  0.83333333333333333,  0.333333333333333333]
    local c2 = SA[ 0.333333333333333333,  0.83333333333333333, -0.166666666666666667]
    
    local U0 = c0[1]*U[i+3] + c0[2]*U[i+2] + c0[3]*U[i+1]
    local U1 = c1[1]*U[i+2] + c1[2]*U[i+1] + c1[3]*U[i]
    local U2 = c2[1]*U[i+1] + c2[2]*U[i]   + c2[3]*U[i-1]


    ## Final computation
    ω0*U0 + ω1*U1 + ω2*U2
end