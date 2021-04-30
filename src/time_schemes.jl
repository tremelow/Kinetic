function swap_pbs!(pbs :: Vector{RelaxPb}, nb_stages)
    if length(pbs) != nb_stages+1
        to_add = (nb_stages + 1) - length(pbs)
        push!(pbs, [copy(pbs[1]) for _ in 1:to_add]... )
    end
    pbs[1], pbs[2:end] = pbs[end], pbs[1:end-1]
end

function imex_bdf1!(pbs :: Vector{RelaxPb}, dt, dx, ε)
    swap_pbs!(pbs,1)

    local Θ = 1.0 / (ε^2 + dt)
    local dx⁻¹ = 1.0 / dx

    ## Compute vⁿ⁺¹
    compute_dxP!(pbs[2], dx⁻¹)
    @. pbs[1].V.data = (1.0 - dt * Θ) * pbs[2].V.data + 
                    dt * Θ * (pbs[2].f(pbs[2].U.data) - pbs[2].dxP )

    ## Compute uⁿ⁺¹
    compute_dxV!(pbs[1], dx⁻¹)
    @. pbs[1].U.data = pbs[2].U.data - dt * pbs[1].dxV

    nothing
end


function imex_bdf2!(pbs :: Vector{RelaxPb}, dt, dx, ε)
    if length(pbs) == 1
        push!(pbs, copy(pbs[1]), copy(pbs[1]))
        imex_bdf1!(pbs, dt, dx, ε)
    else
        swap_pbs!(pbs,2)
        local a0, a1 = 1.333333333333333333, -0.333333333333333333
        local b0, b1 = 1.333333333333333333, -0.666666666666666667
        local c₋₁ = 0.666666666666666667

        local Θ = 1.0 / (ε^2 + c₋₁ * dt)
        local dx⁻¹ = 1.0 / dx

        compute_dxP!(pbs[2], dx⁻¹) # assume pbs[3].dxP is correct
        local U1, U2, U3 = pbs[1].U.data, pbs[2].U.data, pbs[3].U.data
        local V1, V2, V3 = pbs[1].V.data, pbs[2].V.data, pbs[3].V.data
        @. V1 = (1.0 - c₋₁ * dt * Θ) * (a0*V2 + a1*V3) + 
                dt * Θ * (b0 * (pbs[2].f(U2) - pbs[2].dxP ) +
                            b1 * (pbs[3].f(U3) - pbs[3].dxP ) )

        compute_dxV!(pbs[1], dx⁻¹)
        @. U1 = a0*U2 + a1*U3 - c₋₁ * dt * pbs[1].dxV
    end

    nothing
end



function imex_bdf4!(pbs :: Vector{RelaxPb}, dt, dx, ε)
    if length(pbs) == 1
        imex_bdf1!(pbs, dt, dx, ε)
    elseif length(pbs) == 2
        imex_bdf2!(pbs, dt, dx, ε)
    else
        swap_pbs!(pbs,3)
        local a0, b0 =  1.636363636363636364,  1.636363636363636364
        local a1, b1 = -0.818181818181818182, -1.636363636363636364
        local a2, b2 =  0.181818181818181818,  0.545454545454545455
        local c₋₁ = 0.545454545454545455

        local Θ = 1.0 / (ε^2 + c₋₁ * dt)
        local dx⁻¹ = 1.0 / dx

        compute_dxP!(pbs[2], dx⁻¹) # assume pbs[3+].dxP is correct
        local U1, V1 = pbs[1].U.data, pbs[1].V.data
        local U2, V2 = pbs[2].U.data, pbs[2].V.data
        local U3, V3 = pbs[3].U.data, pbs[3].V.data
        local U4, V4 = pbs[4].U.data, pbs[4].V.data
        @. V1 = (1.0 - c₋₁ * dt * Θ) * (a0*V2 + a1*V3 + a2*V4) + 
                dt * Θ * (b0 * (pbs[2].f(U2) - pbs[2].dxP ) +
                            b1 * (pbs[3].f(U3) - pbs[3].dxP ) +
                            b2 * (pbs[4].f(U4) - pbs[4].dxP))

        compute_dxV!(pbs[1], dx⁻¹)
        @. U1 = a0*U2 + a1*U3 + a2*U4 - c₋₁ * dt * pbs[1].dxV
    end

    nothing
end


function imex_bdf3!(pbs :: Vector{RelaxPb}, dt, dx, ε)
    if length(pbs) == 1
        imex_bdf1!(pbs, dt, dx, ε)
    elseif length(pbs) == 2
        imex_bdf2!(pbs, dt, dx, ε)
    elseif length(pbs) == 3
        imex_bdf3!(pbs, dt, dx, ε)
    else
        swap_pbs!(pbs,4)
        local a0, b0 =  1.92,  1.92
        local a1, b1 = -1.44, -2.88
        local a2, b2 =  0.64,  1.92
        local a3, b3 =  0.12, -0.48
        local c₋₁ = 0.48

        local Θ = 1.0 / (ε^2 + c₋₁ * dt)
        local dx⁻¹ = 1.0 / dx

        compute_dxP!(pbs[2], dx⁻¹) # assume pbs[3+].dxP is correct
        local U1, V1 = pbs[1].U.data, pbs[1].V.data
        local U2, V2 = pbs[2].U.data, pbs[2].V.data
        local U3, V3 = pbs[3].U.data, pbs[3].V.data
        local U4, V4 = pbs[4].U.data, pbs[4].V.data
        local U5, V5 = pbs[5].U.data, pbs[5].V.data
        @. V1 = (1.0 - c₋₁ * dt * Θ) * (a0*V2 + a1*V3 + a2*V4 + a3*V5) + 
                dt * Θ * (b0 * (pbs[2].f(U2) - pbs[2].dxP ) +
                            b1 * (pbs[3].f(U3) - pbs[3].dxP ) +
                            b2 * (pbs[4].f(U4) - pbs[4].dxP ) +
                            b3 * (pbs[5].f(U5) - pbs[5].dxP))

        compute_dxV!(pbs[1], dx⁻¹)
        @. U1 = a0*U2 + a1*U3 + a2*U4 - c₋₁ * dt * pbs[1].dxV
    end

    nothing
end