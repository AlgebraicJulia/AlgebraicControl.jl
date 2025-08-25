using AlgebraicControl
using LinearAlgebra
using Convex

# Set up each agent's dynamics: x' = Ax + Bu
dt = 0.1  # Discretization step size
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]
B = [0 0; dt 0; 0 0; 0 dt]
Q = Matrix{Float64}(I(4))
R = Matrix{Float64}(I(2))

sys = LinearSystem(A, B)

# Stage cost function: minimize control effort and deviation from origin
stage_cost(u, x) = quadform(x, Q) + quadform(u, R)
stage_constraints = Function[]
# Create a multi-stage program for a horizon of N steps
N = 100
mpc_program = multi_stage_program(stage_cost, stage_constraints, sys, N)

# Create a proxable MPC program
proxable_mpc = ProxableMPCProgram(mpc_program)

set_x0!(proxable_mpc, [5.0; 0.0; 5.0; 0.0])  # Initial state

xf = prox(proxable_mpc, zeros(4))  # Proximal step towards target state [10; 0; 10; 0]

