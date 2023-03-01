module Optimizers

export optimize, OpenOptimizer, CompositeOptimizer, forward, backward!

using ..Problems
using ..LensCats
using Zygote

mutable struct OpenOptimizer
    dom::Int
    codom::Int
    state::Tuple # (prev_x0, prev_λ0, λ)
    p::Function
end
function OpenOptimizer(p::EqConstrainedProb)
    cd = p.f.dom
    d = p.h.codom
    state = (zeros(cd), zeros(d), zeros(cd))
    return OpenOptimizer(d, cd, state, Problems.open(p))
end

mutable struct CompositeOptimizer
    os::Vector{OpenOptimizer}
end

function forward(o::OpenOptimizer, u)
    return optimize(o.p(o.state[3],u), o.state[1]; λ₀=o.state[2])[1]
end

function backward!(o::OpenOptimizer, u, λ)
    o.state[3] = λ
    x,λ_pb = optimize(o.p(o.state[3], u), o.state[1]; λ₀=o.state[2])
    o.state[1] = x
    o.state[2] = λ_pb
    return λ_pb
end

function forward(os::CompositeOptimizer, u)
    reduce((x,o) -> forward(o, x), os; init=u)
end

function backward!(os::CompositeOptimizer, u, λ)
    # Make a list of each get output
    xs = [u]
    for o in os
        x = xs[end]
        push!(xs, forward(o, x))
    end

    for (o,x) in zip(reverse(os), reverse(xs))
        λ = backward!(o, x, λ)
    end
    return λ
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

function GD(f::Function, x₀::Vector{Float64}; γ=0.001, num_iters=10000, stopping_criterion=0.0001)
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

function ec_uzawa(f, h, x₀, num_constraints; step_size=0.001, max_iters=10000, stopping_criterion=0.0001, λ₀=nothing, tol=0.00001)
    λ = zeros(num_constraints)
    if λ₀ !== nothing
        λ = λ₀
    end
    x = x₀
    γ = step_size
    ϵ = stopping_criterion
    m = num_constraints
    for i in 1:max_iters
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

function uzawa(f,g, x₀, num_constraints; step_size=0.001, max_iters=10000, stopping_criterion=0.0001, μ₀=nothing)
    μ = zeros(num_constraints)
    if μ₀ !== nothing
        μ = μ₀
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