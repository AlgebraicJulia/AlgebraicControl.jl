using AlgebraicControl.Optimizers
using AlgebraicControl.Problems
using AlgebraicControl.BasicCats
using Test
using LinearAlgebra

# Test optimization of an atomic problem
vals = -5.0:5.0
A = rand(2, 5)
B = rand(5, 5)
b = rand(vals, 2)
Q = rand(vals, 5, 5)
#Q = Q'*Q
Q = I(5)
r = rand(vals, 5)

f(x) = x'*Q*x + r'*x
g(x) = A*x - b

p = IneqConstrainedProb(SmoothFunction(5,5,f), SmoothFunction(5,5,g))
o = OpenOptimizer(p)

x_tst1 = forward(o, zeros(5); max_iters=10000)
backward!(o, zeros(5), zeros(5); max_iters=10000)
x_tst2 = primal_value(o)
x_true = optimize(p, zeros(5); max_iters=10000)[1]

@test x_tst1 ≈ x_true
@test x_tst2 ≈ x_true

# Test a composite optimization
f_sub(x) = x'*Q*x + r'*x
f(x) = f_sub(x[1:5]) + f_sub(x[6:end])
h1(x) = A*x-b
h2(x, y) = B*y - x
h(x) = vcat(h1(x[1:5]), h2(x[1:5],x[6:end]))
true_p = EqConstrainedProb(SmoothFunction(10, 1, f), SmoothFunction(10,10,h))
println("Running Composite Optimization...")
true_x, true_μ = optimize(true_p, zeros(10); max_iters=50000)

p1 = EqConstrainedProb(SmoothFunction(5,1,f_sub), SmoothFunction(5,5,h1))
p2 = EqConstrainedProb(SmoothFunction(5,1,f_sub), SmoothFunction(5,5,y->B*y))

println("Running Decomposed Optimization...")
o = compose(OpenOptimizer(p1), OpenOptimizer(p2))
for i in 1:5
    backward!(o, zeros(5), zeros(5))
end

tst_x = vcat(primal_value(o.os[1]), primal_value(o.os[2]))
tst_μ = vcat(dual_value(o.os[1]), dual_value(o.os[2]))

println(norm(tst_x - true_x))
