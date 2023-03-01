using AlgebraicControl.Optimizers
using AlgebraicControl.Problems
using AlgebraicControl.BasicCats
using Test

# Test optimization of an atomic problem
vals = -5.0:5.0
A = rand(vals, 3, 3)
b = rand(vals, 3)
Q = rand(vals, 3, 3)
Q = Q'*Q
r = rand(vals, 3)

f(x) = x'*Q*x + r'*x
h(x) = A*x + b

p = EqConstrainedProb(SmoothFunction(3,3,f), SmoothFunction(3,3,h))
o = OpenOptimizer(p)

x_tst = forward(o, zeros(3))
x_true = optimize(p, zeros(3))[1]

@test x_tst ≈ x_true

