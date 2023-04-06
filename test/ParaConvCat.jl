using AlgebraicControl.ParaConvCat
using Test
using Convex
using LinearAlgebra
using AlgebraicControl.Categories
using SCS

r1 = -5:0.1:5
r2 = 0:5
Q = rand(r1,10,10)
#Q = Q'*Q
#Q = 5*I(10)
Q = diagm(repeat([5],10))
R = rand(r1,8,8)
#R = R'*R
#R = 3*I(8)
R = diagm(repeat([3],8))
A = rand(10,10)
B = rand(10,8)
dim_x = size(Q)[1]
dim_u = size(R)[1]
x0 = repeat([5], 10)

F_impl = (u1, x1, x2) -> ConvexBifunction(
    quadform(x1, Q; assume_psd=true) + quadform(u1[1], R; assume_psd=true),
    [x2 == A*x1 + B*u1[1]]
)
F = OpenParaConvexBifunction(10,10,[8],F_impl)
u1 = Variable(8)
x1 = Variable(10)
x2 = Variable(10)
FCB = F([u1], x1, x2)
prob = to_problem(FCB)
solve!(prob, SCS.Optimizer; silent_solver=true)
o1 = prob.optval

# Test unitality
F_id = compose(ParaConv(), F, id(ParaConv(), 10))
u1 = Variable(8)
x1 = Variable(10)
x2 = Variable(10)
FCB = F_id([u1], x1, x2)
prob = to_problem(FCB)
solve!(prob, SCS.Optimizer; silent_solver=true)
o2 = prob.optval
@test o1 ≈ o2

id_F = compose(ParaConv(), id(ParaConv(), 10), F)
u1 = Variable(8)
x1 = Variable(10)
x2 = Variable(10)
FCB = id_F([u1], x1, x2)
prob = to_problem(FCB)
solve!(prob, SCS.Optimizer; silent_solver=true)
o2 = prob.optval
@test o1 ≈ o2




FF = compose(ParaConv(), F, F)

u1 = Variable(8)
x1 = Variable(10)
x2 = Variable(10)
u2 = Variable(8)
x3 = Variable(10)
FFCB = FF([u1,u2], x1, x3)

prob = to_problem(FFCB)
fix!(x1, repeat([5],10))
fix!(x3, zeros(10))

solve!(prob, SCS.Optimizer; silent_solver=true)
o1 = prob.optval

x1_val = evaluate(x1)
u1_val = evaluate(u1)
x2_val = A*x1_val + B*u1_val
u2_val = evaluate(u2)
x3_val = A*x2_val + B*u2_val

tu1 = Variable(8)
tx1 = Variable(10)
tx2 = Variable(10)
tu2 = Variable(8)
tx3 = Variable(10)
true_prob = minimize(
    quadform(tx1, Q; assume_psd=true) + quadform(tu1, R; assume_psd=true) +
    quadform(tx2, Q; assume_psd=true) + quadform(tu2, R; assume_psd=true),
    [tx2 == A*tx1 + B*tu1, tx3 == A*tx2 + B*tu2]
)
fix!(tx1, repeat([5],10))
fix!(tx3, zeros(10))

solve!(true_prob, SCS.Optimizer, silent_solver=true)
o2 = true_prob.optval

tx1_val = evaluate(tx1)
tu1_val = evaluate(tu1)
tx2_val = A*tx1_val + B*tu1_val
tu2_val = evaluate(tu2)
tx3_val = A*tx2_val + B*tu2_val

@test o1 == o2
@test x1_val == tx1_val
@test x2_val == tx2_val
@test x3_val == tx3_val
@test u1_val == tu1_val
@test u2_val == tu2_val





