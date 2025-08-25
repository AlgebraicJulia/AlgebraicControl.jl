module ParaProbs

export ParaProb, ParaECProb, ParaFunction, param_space, open_prob

using ..Problems
import ..Categories: dom, codom

abstract type ParaProb <: Problem end

struct ParaECProb <: ParaProb
    dom::Int
    codom::Int
    param::Int
    f::Function
    h::Function
end

struct ParaFunction
    dom::Int
    codom::Int
    para::Int
    impl::Function # para × dom → codom
end
(f::ParaFunction)(u::Vector{Float64}, x::Vector{Float64}) = f.impl(u,x)

dom(f::ParaFunction) = f.dom
codom(f::ParaFunction) = f.codom
param_space(f::ParaFunction) = f.para

function open_prob(f::ParaFunction, h::ParaFunction)
    return (λ::Vector{Float64},y::Vector{Float64}) -> begin
        dom(f) == dom(h) || error("Objective and constraint functions must have same domain")
        codom(f) == 1 || error("Objective function must have scalar codomain")
        param_space(f) == param_space(h) || error("Objective and constraint functions must have same parameter space")
        length(λ) == dom(f) || error("Dimension mismatch in open problem construction")
        length(y) == codom(h) || error("Dimension mismatch in open problem construction")
        n_u = param_space(f)
        obj(x) = f.impl(x[1:n_u],x[n_u+1:end]) - λ'*x[n_u+1:end]
        cons(x) = h.impl(x[1:n_u],x[n_u+1:end]) - y
        return ParaECProb(codom(h), dom(f), param_space(f), obj, cons)
    end
end

end

