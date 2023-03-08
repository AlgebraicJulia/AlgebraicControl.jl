using Test
using LinearAlgebra
using AlgebraicControl.Optimizers
using AlgebraicControl.ParaProbs
using AlgebraicControl.ParaOptimizers

A = rand(1:5,3,3)
B = rand(1:3,3,2)
Q = 3*I(3)
R = 5*I(2)

f = ParaFunction(3,1,2,(u,x) -> x'*Q*x + u'*R*u)
h = ParaFunction(3,3,2,(u,x) -> A*x + B*u)
o1 = ParaOptimizer(f, h)
o2 = ParaOptimizer(f, h)

#x = forward(o, ones(3))
#λ = backward!(o, ones(3), zeros(3))

os = compose(o1,o2)

backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
backward!(os, ones(3), zeros(3))
