using Test
using LinearAlgebra
using AlgebraicControl.Optimizers
using AlgebraicControl.LQR
using AlgebraicControl.Problems
using AlgebraicControl.BasicCats

A = [0 1; .01 0]
B = [0.0; 1.0]
Q = [1.0 0; 0 1.0]
R = [1.0]
x₀ = [3.0, 0]

p = generating_problem(A, B, Q, R)

res, λ = optimize(p, zeros(3))

@test length(λ) == 2
@test length(res) == 3

o = OpenOptimizer(p, para=1)
backward!(o, ones(2), zeros(2))

os = compose(o, o)

backward!(os, ones(2), zeros(2))

f(u) = u'*3*u
h(u) = A*x₀ + B*u[1]

#optimize(EqConstrainedProb(SmoothFunction(1,1,f), SmoothFunction(1,2,h)),zeros(1))

ip = init_prob(A, B, Q, R, x₀)
o = OpenOptimizer(ip,para=1)
backward!(o, zeros(2), zeros(0))

LQP = compose(os, o)
backward!(LQP, zeros(2), zeros(0))

#o = LQR_optimizer(A, B, Q, R, x₀, 3)
#optimize!(o)
