module AlgebraicControl

include("ConvexPrograms.jl")

export ConvexProgram, OpenConvexProgram, LinearSystem, single_stage_program, multi_stage_program,
    ProxableMPCProgram, ProxableMPCProgram, set_x0!, prox, prox!

using Convex
using SCS
using LinearAlgebra

import ProximalOperators: prox, prox!

struct LinearSystem
    A::Matrix{Float64}
    B::Matrix{Float64}
    LinearSystem(A::Matrix{Float64}, B::Matrix{Float64}) =
        size(A, 1) == size(B, 1) ? new(A, B) : error("Inconsistent dimensions between A and B matrices.")
end

state_dim(sys::LinearSystem) = size(sys.A, 1)
control_dim(sys::LinearSystem) = size(sys.B, 2)


function (sys::LinearSystem)(x, u)
    return sys.A * x + sys.B * u
end

function single_stage_program(stage_cost::Function, stage_constraints::Vector{Function}, sys::LinearSystem)::OpenConvexProgram
    n = state_dim(sys)
    m = control_dim(sys)
    impl = (u1, x1, x2) -> ConvexProgram(
        stage_cost(u1[1], x1),
        vcat([c(u1[1], x1) for c in stage_constraints]..., x2 == sys(x1, u1[1]))
    )
    return OpenConvexProgram(n, n, [m], impl)
end

function multi_stage_program(stage_cost::Function, stage_constraints::Vector{Function}, sys::LinearSystem, n::Int)::OpenConvexProgram
    one_step = single_stage_program(stage_cost, stage_constraints, sys)
    return compose(one_step, n)
end

struct ProxableMPCProgram
    control_vars::Vector{Variable}
    input_var::Variable
    output_var::Variable
    program
end

ProxableMPCProgram(F::OpenConvexProgram) = begin
    # Make variables of the right dimensions
    x0 = Variable(dom(F))
    xf = Variable(codom(F))
    #y = Variable(codom(F))
    us = [Variable(d) for d in F.internal_vars]
    #γ = Variable(1)
    mpc_program = F(us, x0, xf)
    # Once Convex.jl is fixed, we can make γ and y real variables and fix! them
    # in the prox function instead of rebuilding the program each time. But for now,
    # we suffer....
    proxable_program(y, γ) = minimize( # This is stupid but Convex.jl is broken...
        mpc_program.objective + (1 / (2 * γ)) * LinearAlgebra.dot(xf - y, xf - y),
        mpc_program.constraints
    )
    return ProxableMPCProgram(us, x0, xf, proxable_program)
end

function set_x0!(P::ProxableMPCProgram, x0_val::Vector{Float64})
    fix!(P.input_var, x0_val)
end

function prox(P::ProxableMPCProgram, x, γ=1.0)
    #fix!(P.y, x)
    #fix!(P.γ, γ)
    prob = P.program(x, γ)
    solve!(prob, SCS.Optimizer; silent=true)
    return evaluate(P.output_var), prob.optval
end

function prox!(y, P::ProxableMPCProgram, x, γ=1.0)
    #fix!(P.y, x)
    #fix!(P.γ[], γ)
    prob = P.program(x, γ)
    solve!(prob, SCS.Optimizer; silent=true)
    copy!(y, evaluate(P.output_var))
    return prob.optval
end


end