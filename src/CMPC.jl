using JuMP
using Ipopt
using LinearAlgebra

r1 = -5:0.1:5
r2 = 0:5
Q = rand(r1,10,10)
#Q = Q'*Q
Q = 5*I(10)
R = rand(r1,8,8)
#R = R'*R
R = 3*I(8)
A = rand(10,10)
B = rand(10,8)
dim_x = size(Q)[1]
dim_u = size(R)[1]
x0 = repeat([5], 10)

tp = Model(Ipopt.Optimizer)
@variable(tp, tx1[1:dim_x] >= 0)
@variable(tp, tu1[1:dim_u])
@variable(tp, tx2[1:dim_x] >= 0)
@variable(tp, tu2[1:dim_u])
@objective(tp, Min, tx1'*Q*tx1 + tx2'*Q*tx2 + tu1'*R*tu1 + tu2'*R*tu2)
@constraint(tp, tx1 .== A*tx2 .+ B*tu2)
@constraint(tp, x0 .== A*tx1 .+ B*tu1)
set_optimizer_attribute(tp, "print_level", 0)
#println(tp)
@time optimize!(tp)

mp = Model(Ipopt.Optimizer)
@variable(mp, x1[1:dim_x])
@variable(mp, u1[1:dim_u])
@variable(mp, λ[1:dim_x])
@objective(mp, Min, x1'*Q*x1 + u1'*R*u1 - λ'*x1)
@constraint(mp, x0 .== A*x1 .+ B*u1)
set_optimizer_attribute(mp, "print_level", 0)

sp = Model(Ipopt.Optimizer)
@variable(sp, x2[1:dim_x])
@variable(sp, u2[1:dim_u])
@variable(sp, x1_copy[1:dim_x])
@objective(sp, Min, x2'*Q*x2 + u2'*R*u2)
c = @constraint(sp, x1_copy .== A*x2 .+ B*u2)

set_optimizer_attribute(sp, "print_level", 0)

function solve_masterproblem(model, λ)
    fix.(model[:λ],λ)
    optimize!(model)
    #@assert termination_status(model) == OPTIMAL
    return value.(model[:x1])
end

function solve_subproblem(model, x1)
    fix.(model[:x1_copy],x1)
    optimize!(model)
    #@assert termination_status(model) == OPTIMAL
    return dual.(c)
end

function comp_solve(mp, sp)
    λ = zeros(dim_x)
    for i in 1:100
        x1 = solve_masterproblem(mp, λ)
        λ_new = solve_subproblem(sp, x1)
        if λ_new ≈ λ
            println("Breaking after $i iterations.")
            break
        end
        λ = λ_new
    end
end
@time comp_solve(mp,sp)

x1_true = value.(tp[:tx1])
x2_true = value.(tp[:tx2])
x1 = value.(mp[:x1])
x2 = value.(sp[:x2])
