using Optim
using LinearAlgebra

# Test out constrained optimization with Optim.jl

A = rand(1:5,5,5)
B = rand(1:3,5,3)
Q = 3*I(5)
R = 5*I(3)
x1 = 5
x2 = 6
u1 = 11
u2 = 14

f(in) = in[1:x]'*Q*in[1:x] + 
    in[u:end]'*R*in[u:end]
h!(c, in) = (c .= A*in[1:x] + B*in[u:end];c)
h(in) = A*in[1:x] + B*in[u:end]

obj = TwiceDifferentiable(f, ones(8), autodiff=:forward)
cons = TwiceDifferentiableConstraints(h!, repeat([-Inf],8),repeat([Inf],8),ones(5),ones(5),:forward)

@time sol = optimize(obj, cons, ones(8), IPNewton())
