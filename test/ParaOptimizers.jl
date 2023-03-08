using Test
using LinearAlgebra
using AlgebraicControl.Optimizers
using AlgebraicControl.ParaProbs
using AlgebraicControl.ParaOptimizers

A = rand(1:5,5,5)
B = rand(1:3,5,3)
Q = 3*I(5)
R = 5*I(3)

f = ParaFunction(5,1,3,(u,x) -> x'*Q*x + u'*R*u)
h = ParaFunction(5,5,3,(u,x) -> A*x + B*u)

#=os = ParaOptimizer[]
for i in 1:5
    push!(os, ParaOptimizer(f, h))
end
o = CompositeParaOptimizer(os)

@time optimize!(o, ones(5), zeros(5), num_iters=5)

us = []
xs = []
for o in o.os
    push!(us, param_value(o))
    push!(xs, primal_value(o))
end=#


o1 = ParaOptimizer(f, h)
o2 = ParaOptimizer(f, h)
o3 = ParaOptimizer(f, h)

os = compose(compose(o1,o2), o3)

@time optimize!(os, ones(5), zeros(5), num_iters=15)

u1 = param_value(os.os[1])
u2 = param_value(os.os[2])
u3 = param_value(os.os[3])
x1 = primal_value(os.os[1])
x2 = primal_value(os.os[2])
x3 = primal_value(os.os[3])

res = vcat(u1,u2,u3,x1,x2,x3)

total_f(x) = f(x[1:3], x[10:14]) + f(x[4:6], x[15:19]) + f(x[7:9], x[20:24])
total_h(x) = vcat(h(x[1:3],x[10:14]) - ones(5),
                  h(x[4:6],x[15:19]) - x[10:14],
                  h(x[7:9],x[20:24]) - x[15:19])

@time true_res, _ = ec_uzawa(total_f, total_h, zeros(24), 15; max_iters=100000)
true_u1 = true_res[1:3]
true_u2 = true_res[4:6]
true_u3 = true_res[7:9]
true_x1 = true_res[10:14]
true_x2 = true_res[15:19]
true_x3 = true_res[20:24]

#x = forward(o, ones(3))
#λ = backward!(o, ones(3), zeros(3))



