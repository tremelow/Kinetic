module StiffKinetic
    using Revise

    include("mesh.jl")
    include("problems.jl")
    include("space_derivatives.jl")

    export save_mesh
    export HyperbolicProblem
    export posUpwind!, posWENO3!, posWENO5!
    export negUpwind!, negWENO3!, negWENO5!
    export posUpwind, posWENO3, posWENO5
    export negUpwind, negWENO3, negWENO5
end # module
