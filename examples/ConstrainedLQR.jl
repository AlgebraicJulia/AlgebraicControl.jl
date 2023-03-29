using AlgebraicControl.CMPC
using AlgebraicControl.ParaConvCat
using SCS
using Convex

A = [0 1; .01 0]
B = [0.0; 1.0]
Q = [1.0 0; 0 1.0]
R = 1.0
x₀ = [3.0, 0]
N = 10

cost(u,x) = quadform(x,Q) + square(u)

one_step = one_step_bifunction(cost, A, B)

MPC_bifunc = MPC_bifunction(one_step, N)

us = repeat([Variable(1)], N)
x1 = Variable(2)
x_N = Variable(2)

MPC_prob = to_cvx(MPC_bifunc, us, x1, x_N)

fix!(x1, repeat([5], 2))
solve!(MPC_prob, SCS.Optimizer)
o = MPC_prob.optval


