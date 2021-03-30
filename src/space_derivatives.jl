using StaticArrays
using OffsetArrays
import Base.Iterators: drop, take

include("space_derivatives/upwind.jl")
include("space_derivatives/weno3.jl")
include("space_derivatives/weno5.jl")
