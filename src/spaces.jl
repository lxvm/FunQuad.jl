# integrate smooth spaces using QuadGK (or Clenshaw-Curtis?)
# integrate Jacobi spaces using QuadJGK
# integrate Fourier spaces using AutoSymPTR
# integrate Laurent spaces using MeroQuadGK
# integrate mixed smooth+powerlaw*smooth spaces with generalized Chebyshev
# integrate mixed smooth+logarithm*smooth spaces with generalized Gaussian

# we should cache and reuse the quadrature nodes/weights wherever possible.
# for h-adaptive we should use embedded rules when possible, otherwise h-adaptive estimator
# how to define order of a mixed space?

# function types
struct GaussKronrod end

struct GaussLegendre end

struct LogAtMin end

struct LogAtMax end

# TODO: implement changes of variables with this type
struct Transformed{T,F1,F2,F3}
    fun::T
    x2y::F1
    y2x::F2
    dxdy::F3
end

# quadrature rule traits

struct Embedded end
struct Refined end

EstimatorStyle(::GaussKronrod) = Embedded()
EstimatorStyle(::GaussLegendre) = Refined()
EstimatorStyle(::LogAtMin) = Refined()
EstimatorStyle(::LogAtMax) = Refined()
EstimatorStyle(f::Transformed) = EstimatorStyle(f.fun)

# space constructor
struct FunSpace{T}
    fun
    a::T
    b::T
end

FunSpace(fun, a, b) = FunSpace(fun, promote(a, b)...)
FunSpace(a, b) = FunSpace(GaussKronrod(), a, b)

# default space subdivision methods

split(sp::FunSpace) = _split(sp.fun, sp.a, sp.b)
function _split(sp::GaussKronrod, a, b)
    mid = (a+b)/2
    sp1 = FunSpace(sp, oftype(mid, a), mid)
    sp2 = FunSpace(sp, mid, oftype(mid, b))
    return sp1, sp2
end
# TODO: make GaussLegendre do lazy evaluation of Kronrod points (i.e. embedded)
function _split(sp::GaussLegendre, a, b)
    mid = (a+b)/2
    sp1 = FunSpace(sp, oftype(mid, a), mid)
    sp2 = FunSpace(sp, mid, oftype(mid, b))
    return sp1, sp2
end

function _split(sp::LogAtMin, a, b)
    mid = (a+b)/2
    sp1 = FunSpace(sp, oftype(mid, a), mid)
    sp2 = FunSpace(GaussLegendre(), mid, oftype(mid, b))
    return sp1, sp2
end
function _split(sp::LogAtMax, a, b)
    mid = (a+b)/2
    sp1 = FunSpace(GaussLegendre(), oftype(mid, a), mid)
    sp2 = FunSpace(sp, mid, oftype(mid, b))
    return sp1, sp2
end

function cachednodes(sp::FunSpace{T}, n::Integer)::Vector{real(float(T))} where {T}
    # intentional dynamic dispatch, since type check will guarantee stability
    return first(_cachedrule(sp.fun, T, Int(n)))
end

function _cachedrule(::GaussKronrod, ::Type{T}, n::Int) where {T}
    return QuadGK.cachedrule(T, n)
end

const gaussrulecache = Dict{Type,Dict}()

@generated function _cachedgaussrule(::Type{T}, n::Int) where {T}
    cache = haskey(gaussrulecache, T) ? gaussrulecache[T] : (gaussrulecache[T] = Dict{Int,NTuple{2,Vector{T}}}())
    :(haskey($cache, n) ? $cache[n] : ($cache[n] = QuadGK.gauss(T, n)))
end

function _cachedrule(::GaussLegendre, ::Type{T}, n::Int) where {T}
    return _cachedgaussrule(T, n)
end

const logrulecache = Dict{Type,Dict}()

function logrule(::Type{T}, n::Int) where {T}
    try
        x, w = GeneralizedGaussianQuadrature.generalizedquadrature(n)
        return Vector{T}(x), Vector{T}(w)
    catch
        throw(ArgumentError("order $n quadrature not available, only 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25 and 30"))
    end
end

@generated function _cachedlogrule(::Type{T}, n::Int) where {T}
    cache = haskey(logrulecache, T) ? logrulecache[T] : (logrulecache[T] = Dict{Int,NTuple{2,Vector{T}}}())
    :(haskey($cache, n) ? $cache[n] : ($cache[n] = logrule(T, n)))
end

function _cachedrule(::LogAtMin, ::Type{T}, n::Int) where {T}
    return _cachedlogrule(typeof(real(float(one(T)))), n)
end

function _cachedrule(::LogAtMax, ::Type{T}, n::Int) where {T}
    return _cachedrule(LogAtMin(), T, n)
end

struct Panel{TX,TI,TE}
    space::FunSpace{TX}
    I::TI
    E::TE
    children
end

Base.isless(p1::Panel, p2::Panel) = isless(p1.E, p2.E)

function _evalrule(f::F, ::Embedded, sp::GaussKronrod, a, b, n, nrm) where {F}
    x, w, gw = _cachedrule(sp, typeof(a), n)
    seg = QuadGK.evalrule(f, a, b, x, w, gw, nrm)
    return Panel(FunSpace(sp, seg.a, seg.b), seg.I, seg.E, nothing)
end

function evalrule(f, sp::GaussLegendre, a, b, n)
    s = (b-a)/2
    x, w = _cachedrule(sp, typeof(s), n)
    I = w[1] * f(a + (1 + x[1]) * s)
    for i in 2:n
        I += w[i] * f(a + (1 + x[i]) * s)
    end
    return Panel(FunSpace(sp, oftype(s, a), oftype(s, b)), I*s, nothing, nothing)
end

function evalrule(f, sp::LogAtMin, a, b, n)
    s = float(b-a)
    x, w = _cachedrule(sp, typeof(s), n)
    I = w[1] * f(a + x[1] * s)
    for i in 2:n
        I += w[i] * f(a + x[i] * s)
    end
    return Panel(FunSpace(sp, oftype(s, a), oftype(s, b)), I * s, nothing, nothing)
end
# need to reverse the nodes, weights
function evalrule(f, sp::LogAtMax, a, b, n)
    s = float(b-a)
    x, w = _cachedrule(sp, typeof(s), n)
    I = w[1] * f(a + (1 - x[1]) * s)
    for i in 2:n
        I += w[i] * f(a + (1 - x[i]) * s)
    end
    return Panel(FunSpace(sp, oftype(s, a), oftype(s, b)), I * s, nothing, nothing)
end

function _evalrule(f::F, ::Refined, fun, a_, b_, n, nrm) where {F}
    a, b = float(a_), float(b_)
    sp = FunSpace(fun, a, b)
    p0 = evalrule(f, fun, a, b, n)
    sp1, sp2 = split(sp)
    p1 = evalrule(f, sp1.fun, sp1.a, sp1.b, n)
    p2 = evalrule(f, sp2.fun, sp2.a, sp2.b, n)
    I = p1.I + p2.I
    E = nrm(I - p0.I)
    return Panel(sp, I, E, (p1, p2))
end

function evalrule(f::F, sp::FunSpace, n, nrm) where {F}
    return _evalrule(f, EstimatorStyle(sp.fun), sp.fun, sp.a, sp.b, n, nrm)
end

function refine2(f::F, panel::P, n::Integer, nrm)::Tuple{P,P} where {F,P}
    return refine2_(f, EstimatorStyle(panel.space.fun), panel.space, panel.children, n, nrm)
end

function refine2_(f::F, ::Embedded, sp, ::Nothing, n, nrm) where {F}
    sp1, sp2 = split(sp)
    p1 = evalrule(f, sp1, n, nrm)
    p2 = evalrule(f, sp2, n, nrm)
    return p1, p2
end
function refine2_(f::F, ::Refined, sp, (p1, p2), n, nrm) where {F}
    sp11, sp12 = split(p1.space)
    p11 = evalrule(f, sp11.fun, sp11.a, sp11.b, n)
    p12 = evalrule(f, sp12.fun, sp12.a, sp12.b, n)
    I1 = p11.I + p12.I
    E1 = nrm(p1.I - I1)
    p1_ = Panel(p1.space, I1, E1, (p11, p12))
    sp21, sp22 = split(p2.space)
    p21 = evalrule(f, sp21.fun, sp21.a, sp21.b, n)
    p22 = evalrule(f, sp22.fun, sp22.a, sp22.b, n)
    I2 = p21.I + p22.I
    E2 = nrm(p2.I - I2)
    p2_ = Panel(p2.space, I2, E2, (p21, p22))
    return p1_, p2_
end
