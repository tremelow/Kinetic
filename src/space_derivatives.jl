using StaticArrays
using OffsetArrays
import Base.Iterators: drop, take
import LinearAlgebra: dot

include("space_derivatives/upwind.jl")
include("space_derivatives/weno3.jl")
include("space_derivatives/weno5.jl")
include("space_derivatives/second_order.jl")