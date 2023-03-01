module Optimizers

export optimize

using ..Problems
using ..LensCats
using Zygote


#OpenOptimizer = Lens{Int, Function}
function OpenOptimizer(p::EqConstrainedProb)
    @assert p.f.dom == p.h.dom
    n = p.f.dom
    m = p.h.codom

    
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