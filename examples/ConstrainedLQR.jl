using AlgebraicControl.CMPC
#using AlgebraicControl.ParaConvCat
using SCS
#using Convex

A = [0 1; .01 0]
B = [0; 1]
Q = 5*[1.0 0; 0 1.0]
R = 3.0
x₀ = [3, 1]
N = 10

cost(uₖ,xₖ) = quadform(xₖ,Q) + R*square(uₖ)
set_constraints(uₖ,xₖ) = [
    uₖ <= 1, uₖ >= -1,
    xₖ[1] <= 3, xₖ[1] >= -3,
    xₖ[2] <= 2, xₖ[2] >= -2
]
dynamics(uₖ,xₖ,xₖplus1) = xₖplus1 == A*xₖ + B*uₖ

one_step = one_step_bifunction(cost, set_constraints, A, B)

MPC_bifunc = MPC_bifunction(one_step, N)

us = [Variable(1) for i in 1:N]
#x1 = Variable(2)
#x1 = repeat([1],2)
x_N = Variable(2)
#x_N = zeros(2)

MPC_prob = to_cvx(MPC_bifunc, us, x₀, x_N)

#fix!(x1, repeat([1], 2))
solve!(MPC_prob, SCS.Optimizer)
o = MPC_prob.optval

function simulate(A, B, x₀, us, N)
    x = x₀
    for i in 1:N
        x = A*x + B*us[i]
    end
    return x
end
