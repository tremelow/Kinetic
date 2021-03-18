# @doc raw"""
#     HyperbolicProblem(ε,p,f)

# Defines a structure to represent the problem 
# ```math
# \begin{cases}
#     ∂ₜ u + ∂ₓ v = 0 , \\
#     ∂ₜ v + \frac{1}{ε^{2α}} ∂ₓ p(u) = -\frac{1}{ε^{1+α}} ( v - f(u) ) .
# \end{cases}
# ```
# """
struct HyperbolicProblem{εT <: AbstractFloat}
    ε :: εT
    p :: Function
    f :: Function
    
    function HyperbolicProblem(ε :: εT, p :: Function, 
                               f :: Function) where εT <: AbstractFloat
        new{εT}(ε, p, f)
    end

    HyperbolicProblem() = new{Float64}(1.0, u -> u, u -> 0.0)

    function HyperbolicProblem(ε :: εT) where εT <: AbstractFloat
        new{εT}(ε, u -> u, u -> 0.0)
    end
end
