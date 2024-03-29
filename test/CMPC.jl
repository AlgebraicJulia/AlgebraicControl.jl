using Test
using AlgebraicControl.CMPC
using AlgebraicControl.ParaConvCat
using Convex
using SCS

A = [0 1; .01 0]
B = [0; 1]
Q = 5*[1.0 0; 0 1.0]
R = 3.0
x₀ = [3, 1]
N = 10
dim_x = 2
dim_u = 1
ϵ = 0.0001

cost(uₖ,xₖ) = quadform(xₖ,Q) + R*square(uₖ)
dynamics(uₖ,xₖ) = A*xₖ + B*uₖ
set_constraints(uₖ,xₖ) = [
    uₖ <= 1, uₖ >= -1,
    xₖ[1] <= 3, xₖ[1] >= -3,
    xₖ[2] <= 2, xₖ[2] >= -2
]

one_step = one_step_bifunction(dim_x,dim_u,cost, set_constraints, dynamics)
MPC_bifunc = MPC_bifunction(one_step, N)

us = [Variable(1) for i in 1:N-1]
x_N = Variable(2)

MPC_prob = to_cvx(MPC_bifunc, us, x₀, x_N)

solve!(MPC_prob, SCS.Optimizer; silent_solver=true)

us_sol = evaluate.(us)

function simulate(A, B, x₀, us, N)
    x = x₀
    for i in 1:N
        x = A*x + B*us[i]
    end
    return x
end

# Test that simulate comes close to x_N
res = simulate(A, B, x₀, us_sol, N-1)
@test norm(res - evaluate(x_N)) < ϵ

### Compare to hand written implementation
xs = Variable(2, N)
us = Variable(1, N-1)

constraints = Constraint[
    xs[:, i+1] == A*xs[:,i] + B*us[:,i] for i in 1:N-1
]

for i in 1:N-1
    push!(constraints, xs[:,i][1] <= 3)
    push!(constraints, xs[:,i][1] >= -3)
    push!(constraints, xs[:,i][2] <= 2)
    push!(constraints, xs[:,i][2] >= -2)
    push!(constraints, us[:,i] <= 1)
    push!(constraints, us[:,i] >= -1)
end

push!(constraints, xs[:,1] == x₀)

objective = sum([quadform(xs[:,i], Q) + R*square(us[:,i]) for i in 1:N-1])

prob = minimize(objective, constraints)

solve!(prob, SCS.Optimizer; silent_solver=true)

true_us_sol = evaluate.(us)

# Test that hand written implementation gets same solution
# as automated implementation.
@test norm(true_us_sol - us_sol) < ϵ