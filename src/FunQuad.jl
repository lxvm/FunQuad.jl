module FunQuad

using LinearAlgebra:norm
using DataStructures: heapify!, heappop!, heappush!
# using GeneralizedChebyshevQuadrature: generalizedchebyshev
import GeneralizedGaussianQuadrature
import QuadGK
# using QuadJGK: quadjgk
import Base.Order: Reverse

export GaussKronrod, GaussLegendre, LogAtMin, LogAtMax, FunSpace
export funquadh, funquadh_count

include("spaces.jl")
include("hadapt.jl")
include("padapt.jl")

end # module FunQuad
