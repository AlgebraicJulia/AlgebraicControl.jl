module Problems

export UnconstrainedProb, EqConstrainedProb, IneqConstrainedProb, ConstrainedProb

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

end
