using AlgebraicOptimization
using AlgebraicControl
using Convex
using LinearAlgebra
using ProximalAlgorithms

# Set up each agent's dynamics: x' = Ax + Bu
dt = 0.1  # Discretization step size
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]
B = [0 0; dt 0; 0 0; 0 dt]
C = [1.0 0 0 0; 0 0 1.0 0] # Output agent's position
Q = Matrix{Float64}(I(4))
Q[1, 1] = 0.0
Q[3, 3] = 0.0 # No objective on positions, just make velocities go to 0
R = Matrix{Float64}(I(2))

sys = LinearSystem(A, B)

# Stage cost function: minimize control effort and deviation from origin
stage_cost(u, x) = quadform(x, Q) + quadform(u, R)
stage_constraints = Function[]
# Create a multi-stage program for a horizon of N steps
N = 10
mpc_program = multi_stage_program(stage_cost, stage_constraints, sys, N)

# Create a proxable MPC program
agent1_obj, agent2_obj, agent3_obj = ProxableMPCProgram(mpc_program), ProxableMPCProgram(mpc_program), ProxableMPCProgram(mpc_program)

set_x0!(agent1_obj, rand(4))
set_x0!(agent2_obj, rand(4))
set_x0!(agent3_obj, rand(4))

s = @cellular_sheaf C begin
    x::Stalk{4}, y::Stalk{4}, z::Stalk{4}

    C(x) == C(y)
    C(x) == C(z)
    C(y) == C(z)
end


p = HomologicalProgram([agent1_obj, agent2_obj, agent3_obj], s)

solve(p, ProximalAlgorithms.DouglasRachford(maxit=10))

sim_length = 100

for i in 1:sim_length
    solve(p, ProximalAlgorithms.DouglasRachford(maxit=10))
    x1_curr = evaluate(agent1_obj.input_var)
    x2_curr = evaluate(agent2_obj.input_var)
    x3_curr = evaluate(agent3_obj.input_var)

    u1_curr = evaluate(agent1_obj.control_vars[1])
    u2_curr = evaluate(agent2_obj.control_vars[1])
    u3_curr = evaluate(agent3_obj.control_vars[1])

    x1_next = sys(x1_curr, u1_curr)
    x2_next = sys(x2_curr, u2_curr)
    x3_next = sys(x3_curr, u3_curr)

    set_x0!(agent1_obj, x1_next)
    set_x0!(agent2_obj, x2_next)
    set_x0!(agent3_obj, x3_next)
end

x1_f = evaluate(agent1_obj.input_var)
x2_f = evaluate(agent2_obj.input_var)
x3_f = evaluate(agent3_obj.input_var)



