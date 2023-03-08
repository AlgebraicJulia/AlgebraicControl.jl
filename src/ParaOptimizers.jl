module ParaOptimizers

export ParaOptimizer, param_value, CompositeParaOptimizer

using ..Optimizers
import ..Optimizers: optimize, forward, backward!, compose, primal_value,
    optimize!
using ..ParaProbs
import ..Categories: dom, codom

mutable struct ParaOptimizer
    dom::Int
    codom::Int
    para::Int
    warm_starts::Vector{Vector} # {prev_u0, prev_x0, prev_λ0}
    state::Vector
    p::Function
end
function ParaOptimizer(f::ParaFunction, h::ParaFunction)
    d = codom(h)
    cd = dom(f)
    para = param_space(f)
    warm_starts = [zeros(para), zeros(cd), zeros(d)]
    state = zeros(cd)
    p = open_prob(f, h)
    return ParaOptimizer(d, cd, para, warm_starts, state, p)
end

primal_value(o::ParaOptimizer) = o.warm_starts[2]
dual_value(o::ParaOptimizer) = o.warm_starts[3]
param_value(o::ParaOptimizer) = o.warm_starts[1]

mutable struct CompositeParaOptimizer
    os::Vector{ParaOptimizer}
end

dom(o::CompositeParaOptimizer) = o.os[1].dom
codom(o::CompositeParaOptimizer) = o.os[end].codom

function compose(o1::ParaOptimizer, o2::ParaOptimizer)
    @assert o1.codom == o2.dom
    return CompositeParaOptimizer([o1,o2])
end
function compose(o1::CompositeParaOptimizer, o2::ParaOptimizer)
    @assert o1.os[end].codom == o2.dom
    return CompositeParaOptimizer(vcat(o1.os, o2))
end
function compose(o1::ParaOptimizer, o2::CompositeParaOptimizer)
    @assert o1.codom == o2.os[1].dom
    return CompositeParaOptimizer(vcat(o1, o2.os))
end
function compose(o1::CompositeParaOptimizer, o2::CompositeParaOptimizer)
    @assert o1.os[end].codom == o2.os[1].dom
    return CompositieParaOptimizer(vcat(o1.os, o2.os))
end

function optimize(p::ParaECProb, x₀; kwargs...)
    return ec_uzawa(p.f, p.h, x₀, p.dom; kwargs...)
end

function forward(o::ParaOptimizer, y::Vector)
    length(y) == o.dom || error("Type error, input to forward not domain type")
    prob = o.p(o.state, y)
    return optimize(prob, vcat(o.warm_starts[1], o.warm_starts[2]); dual_init=o.warm_starts[3])[1][o.para+1:end]
end

function backward!(o::ParaOptimizer, y::Vector, λ::Vector)
    length(y) == o.dom || error("Type error, input to backward! not domain type")
    length(λ) == o.codom || error("Type error, input to backward! not param type")
    o.state = λ
    prob = o.p(o.state, y)
    res_x, res_λ = optimize(prob, vcat(o.warm_starts[1], o.warm_starts[2]); dual_init=o.warm_starts[3])
    u = res_x[1:o.para]
    x = res_x[o.para+1:end]
    length(x) == o.codom || error("Type error, result doesn't match codom type")
    o.warm_starts = [u,x,res_λ]
    return res_λ
end

function forward(os::CompositeParaOptimizer, y::Vector)
    length(y) == dom(os) || error("Type error, input to forward not domain type")
    reduce((x,o) -> forward(o,x), os; init=y)
end

function backward!(os::CompositeParaOptimizer, y::Vector, λ::Vector)
    length(y) == dom(os) || error("Type error, input to backward! not domain type")
    length(λ) == codom(os) || error("Type error, input to backward! not param type")
    λ_new = λ
    for i in length(os.os):-1:1
        if i != 1
            λ_new = backward!(os.os[i], primal_value(os.os[i-1]), λ_new)
        else
            λ_new = backward!(os.os[i], y, λ_new)
        end
    end
    return λ_new
end

function optimize!(o::CompositeParaOptimizer, y::Vector, λ::Vector; num_iters=10)
    for i in 1:num_iters
        backward!(o, y, λ)
    end
end

end # module