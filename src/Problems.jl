module Problems

export UnconstrainedProb, EqConstrainedProb, IneqConstrainedProb, ConstrainedProb,
    open

using ..BasicCats

abstract type Problem end

struct UnconstrainedProb <: Problem
    f::SmoothFunction
end

struct EqConstrainedProb <: Problem
    f::SmoothFunction
    h::SmoothFunction
end

struct IneqConstrainedProb <: Problem
    f::SmoothFunction
    g::SmoothFunction
end

struct ConstrainedProb <: Problem
    f::SmoothFunction
    g::SmoothFunction
    h::SmoothFunction
end

# "Open" a problem by making it responsive to inputs
function open(p::EqConstrainedProb) 
    m = p.h.codom
    n = p.f.dom
    return (λ,u) -> EqConstrainedProb(
        SmoothFunction(n, 1, x -> p.f(x) - λ'*x), 
        SmoothFunction(n, m, x -> p.h(x) - u))
end

end
