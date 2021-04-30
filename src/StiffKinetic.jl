module StiffKinetic
    using Revise

    include("problems.jl")
    include("space_derivatives.jl")
    include("time_schemes.jl")

    # export save_mesh
    export RelaxPb
    export posUpwind!, posWENO3!, posWENO5!
    export negUpwind!, negWENO3!, negWENO5!
    export posUpwind, posWENO3, posWENO5
    export negUpwind, negWENO3, negWENO5
    export fd_diff2_ord2!, fd_diff2_ord4!, fd_diff2_ord6!

    export swap_pbs!
    export imex_bdf1!, imex_bdf2!, imex_bdf3!, imex_bdf4!
    
end # module
