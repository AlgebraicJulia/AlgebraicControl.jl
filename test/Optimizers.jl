using AlgebraicControl.Optimizers
using AlgebraicControl.Problems
using AlgebraicControl.BasicCats
using Test
using LinearAlgebra

# Test optimization of an atomic problem
vals = -5.0:5.0
A = rand(vals, 3, 3)
B = rand(vals, 3, 3)
b = rand(vals, 3)
Q = rand(vals, 3, 3)
Q = Q'*Q
r = rand(vals, 3)

f(x) = x'*Q*x + r'*x
g(x) = A*x - b

p = IneqConstrainedProb(SmoothFunction(3,3,f), SmoothFunction(3,3,g))
o = OpenOptimizer(p)

x_tst1 = forward(o, zeros(3); max_iters=10000)
backward!(o, zeros(3), zeros(3); max_iters=10000)
x_tst2 = primal_value(o)
x_true = optimize(p, zeros(3); max_iters=10000)[1]

@test x_tst1 ≈ x_true
@test x_tst2 ≈ x_true

# Test a composite optimization
f_sub(x) = x'*Q*x + r'*x
f(x) = f_sub(x[1:3]) + f_sub(x[4:end])
g1(x) = A*x-b
g2(x, y) = B*y - x
g(x) = vcat(g1(x[1:3]), g2(x[1:3],x[4:end]))
true_p = IneqConstrainedProb(SmoothFunction(6, 1, f), SmoothFunction(6,6,g))
true_x, true_μ = optimize(true_p, zeros(6); max_iters=100000)

p1 = IneqConstrainedProb(SmoothFunction(3,1,f_sub), SmoothFunction(3,3,g1))
p2 = IneqConstrainedProb(SmoothFunction(3,1,f_sub), SmoothFunction(3,3,y->B*y))

o = compose(OpenOptimizer(p1), OpenOptimizer(p2))
for i in 1:10
    backward!(o, zeros(3), zeros(3))
end

tst_x = vcat(primal_value(o.os[1]), primal_value(o.os[2]))
tst_μ = vcat(dual_value(o.os[1]), dual_value(o.os[2]))

println(norm(tst_x - true_x))
