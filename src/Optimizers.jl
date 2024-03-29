module Optimizers

export optimize, optimize!, OpenOptimizer, CompositeOptimizer, forward, backward!,
    primal_value, dual_value, params, compose, primal_values, ec_uzawa

using ..Problems
using ..LensCats
using Zygote

global debug = false

mutable struct OpenOptimizer
    dom::Int
    codom::Int
    state::Vector{Vector} # (prev_x0, prev_λ0, λ)
    para::Int
    p::Function
end
function OpenOptimizer(p::Problem; para=0)
    cd = objective(p).dom-para
    d = constraints(p).codom
    state = [zeros(cd+para), zeros(d), zeros(cd)]
    return OpenOptimizer(d, cd, state, para, Problems.open(p, para=para))
end
primal_value(o::OpenOptimizer) = o.state[1]
dual_value(o::OpenOptimizer) = o.state[2]
params(o::OpenOptimizer) = o.state[3]

mutable struct CompositeOptimizer
    os::Vector{OpenOptimizer}
end
primal_values(os::CompositeOptimizer) = [primal_value(o) for o in os.os]


function compose(o1::OpenOptimizer, o2::OpenOptimizer)
    @assert o1.codom == o2.dom
    return CompositeOptimizer([o1,o2])
end
function compose(o1::CompositeOptimizer, o2::OpenOptimizer)
    @assert o1.os[end].codom == o2.dom
    return CompositeOptimizer(vcat(o1.os, o2))
end
function compose(o1::OpenOptimizer, o2::CompositeOptimizer)
    @assert o1.codom == o2.os[1].dom
    return CompositeOptimizer(vcat(o1, o2.os))
end
function compose(o1::CompositeOptimizer, o2::CompositeOptimizer)
    @assert o1.os[end].codom == o2.os[1].dom
    return CompositieOptimizer(vcat(o1.os, o2.os))
end

function forward(o::OpenOptimizer, u; kwargs...)
    return optimize(o.p(o.state[3],u), o.state[1]; dual_init=o.state[2], kwargs...)[1]
end

function backward!(o::OpenOptimizer, u, λ; kwargs...)
    o.state[3] = λ
    x,λ_pb = optimize(o.p(o.state[3], u), o.state[1]; dual_init=o.state[2], kwargs...)
    o.state[1] = x
    o.state[2] = λ_pb
    return λ_pb
end

function forward(os::CompositeOptimizer, u)
    reduce((x,o) -> forward(o, x), os, init=u)
end

function backward!(os::CompositeOptimizer, u, λ)
    # Make a list of each get output
    xs = vcat(u, primal_values(os))
    #=for o in os.os
        x = xs[end]
        push!(xs, forward(o, x))
    end=#

    for (o,x) in zip(reverse(os.os), reverse(xs))
        λ = backward!(o, x[1:end-o.para], λ)
    end
    return λ
end

function optimize!(os::CompositeOptimizer; iters=10)
    m = os.os[1].dom
    n = os.os[end].codom
    for i in 1:iters
        backward!(os, zeros(m), zeros(n))
    end
end

function optimize(p::UnconstrainedProb, x₀; kwargs...)
    return GD(p.f.impl, x₀; kwargs...)
end

function optimize(p::EqConstrainedProb, x₀; kwargs...)
    m = p.h.codom
    return ec_uzawa(p.f.impl, p.h.impl, x₀, m; kwargs...)
end

function optimize(p::IneqConstrainedProb, x₀; kwargs...)
    m = p.g.codom
    return uzawa(p.f.impl, p.g.impl, x₀, m; kwargs...)
end

# Project a vector onto the non-negative orthant
project(y::Vector{Float64}) = map(x -> x < 0 ? 0 : x, y)

function GD(f::Function, x₀::Vector{Float64}; γ=0.001, max_iters=10000, stopping_criterion=0.0001)
    x = x₀
    ϵ = stopping_criterion
    for k in 1:num_iters
        x_new = x - γ*gradient(f, x)[1]
        if abs(f(x_new) - f(x)) < ϵ
            println("Converged in $k iterations.")
            break
        end
    end
    return x
end

function ec_uzawa(f, h, x₀, num_constraints; step_size=0.001, max_iters=5000, stopping_criterion=0.0001, dual_init=nothing, tol=0.0001)
    λ = zeros(num_constraints)
    if dual_init !== nothing
        λ = dual_init
    end
    x = x₀
    γ = step_size
    ϵ = stopping_criterion
    m = num_constraints
    for i in 1:max_iters
        if debug
            println("Iteration $i")
            println("x = $x")
            println(gradient(f, x)[1])
            println(jacobian(h,x)[1])
        end
        x_new = x - γ*(gradient(f, x)[1] + jacobian(h, x)[1]'*λ)
        λ_new = λ + γ*h(x)
        if abs(f(x_new) - f(x)) < ϵ && reduce(&, repeat([-tol], m) .< h(x) .< repeat([tol], m))
            println("Converged in $i iterations.")
            break
        end
        x = x_new
        λ = λ_new
    end
    return x,λ
end

function uzawa(f,g, x₀, num_constraints; step_size=0.001, max_iters=10000, stopping_criterion=0.0001, dual_init=nothing)
    μ = zeros(num_constraints)
    if dual_init !== nothing
        μ = dual_init
    end
    x = x₀
    γ = step_size
    ϵ = stopping_criterion
    n = length(x)
    m = num_constraints
    for i in 1:max_iters
        x_new = x - γ*(gradient(f, x)[1] + jacobian(g, x)[1]'*μ)
        μ_new = project(μ + γ*g(x))
        if abs(f(x_new) - f(x)) < ϵ && sum(g(x) .<= zeros(m)) == m
            println("Converged in $i iterations.")
            break
        end
        x = x_new
        μ = μ_new
    end
    return x,μ
end


end