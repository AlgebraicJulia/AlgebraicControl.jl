module Problems

export Problem, UnconstrainedProb, EqConstrainedProb, IneqConstrainedProb, ConstrainedProb,
    open, objective, constraints

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

objective(p::EqConstrainedProb) = p.f
objective(p::IneqConstrainedProb) = p.f
objective(p::ConstrainedProb) = p.f

constraints(p::EqConstrainedProb) = p.h
constraints(p::IneqConstrainedProb) = p.g
#TODO: constraints(p::ConstrainedProb)


# "Open" a problem by making it responsive to inputs
function open(p::EqConstrainedProb; para=0) 
    m = p.h.codom
    n = p.f.dom
    n_x = n - para
    return (λ,u) -> EqConstrainedProb(
        SmoothFunction(n, 1, x -> p.f(x) - λ'*x[1:n_x]), 
        SmoothFunction(n, m, x -> p.h(x) - u))
end

function open(p::IneqConstrainedProb) 
    m = p.g.codom
    n = p.f.dom
    return (μ,u) -> IneqConstrainedProb(
        SmoothFunction(n, 1, x -> p.f(x) - μ'*x), 
        SmoothFunction(n, m, x -> p.g(x) - u))
end

end
