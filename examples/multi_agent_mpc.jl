using AlgebraicOptimization
using AlgebraicControl
using Convex
using LinearAlgebra
using ProximalAlgorithms
using MatrixEquations
using SparseArrays
using Plots
using BlockArrays
using ProximalOperators
using ForwardDiff

# Set up each agent's dynamics: x' = Ax + Bu
dt = 0.1  # Discretization step size
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]
B = [0 0; dt 0; 0 0; 0 dt]
C = [1.0 0 0 0; 0 0 1.0 0] # Output agent's position
#C = [0 1.0 0 0; 0 0 0 1.0] # Velocity output

Q = Matrix{Float64}(I(4))
#Q[1, 1] = 0.0
#Q[3, 3] = 0.0
R = Matrix{Float64}(I(2))

Aa = [A zeros(4, 2); C I(2)]
Ba = [B; zeros(2, 2)]
Qa = Matrix{Float64}(I(6))
Qa[2, 2] = 0
Qa[4, 4] = 0
Qa[5, 5] = 100
Qa[6, 6] = 100

X, _ = ared(Aa, Ba, R, Qa)

function feedback_gain(A, B, R, X)
    return -inv(R + B' * X * B) * (B' * X * A)
end

#F = -inv(R + Ba' * X * Ba) * (Ba' * X * Aa)
F = feedback_gain(Aa, Ba, R, X)
function optimal_control(A, B, F, state_dim, x_curr, y_ref)
    u = F * x_curr
    y_ref_full = vcat(zeros(state_dim), y_ref)
    x_next = A * x_curr + B * u - y_ref_full
    return u, x_next
end

function compute_trajectory(A, B, F, x0, n, state_dim; reference_traj=nothing)
    traj = [x0]
    us = Vector{Float64}[]
    for i in 1:n-1
        x_curr = traj[end]
        u = F * x_curr
        y_ref = isnothing(reference_traj) ? zeros(length(x0)) : vcat(zeros(state_dim), reference_traj[i])
        x_next = A * x_curr + B * u - y_ref
        push!(traj, x_next)
        push!(us, u)
    end
    return traj, us
end
reference_traj = vcat(repeat([[1.0, 1.0]], 100), repeat([[-1.0, 1.0]], 100))
N = 200
#reference_traj = [[sin(0.05 * i) * 5.0, cos(0.05 * i) * 5.0] for i in 1:N]
traj, us = compute_trajectory(Aa, Ba, F, vcat(rand(4), zeros(2)), N, 4, reference_traj=reference_traj)

function plot_traj_2d(traj, pos_map, title_str)
    x_pos = [(pos_map*traj[i])[1] for i in eachindex(traj)]
    y_pos = [(pos_map*traj[i])[2] for i in eachindex(traj)]
    plt = plot(x_pos, y_pos, title=title_str, xlabel="x", ylabel="y", legend=false)
    return plt
end

function plot_traj_2d!(plt, traj, pos_map)
    x_pos = [(pos_map*traj[i])[1] for i in eachindex(traj)]
    y_pos = [(pos_map*traj[i])[2] for i in eachindex(traj)]
    plot!(plt, x_pos, y_pos)
end


function animate_trajs_2d(trajs, title_str, labels, xlims, ylims; fps=10)
    anim = @animate for t in 1:length(trajs[1])
        plt = plot(title=title_str, xlabel="x", ylabel="y", xlims=xlims, ylims=ylims, legendposition=:topright)
        for (traj, label) in zip(trajs, labels)
            scatter!(plt, [traj[t][1]], [traj[t][2]], ms=5, label=label)
        end
    end
    return gif(anim, "anim.gif", fps=fps)
end

function animate_trajs_2d(trajs, pos_map, title_str, labels, xlims, ylims; fps=10)
    anim = @animate for t in 1:length(trajs[1])
        plt = plot(title=title_str, xlabel="x", ylabel="y", xlims=xlims, ylims=ylims, legendposition=:topright)
        for (traj, label) in zip(trajs, labels)
            x = pos_map * traj[t]
            scatter!(plt, [x[1]], [x[2]], ms=5, label=label)
        end
    end
    return gif(anim, "anim.gif", fps=fps)
end


Ca = [[1.0 0 0 0; 0 0 1.0 0] zeros(2, 2)]
plot_traj_2d(traj, Ca, "LQ Tracking Trajectory")

#Ca = [[1.0 0 0 0; 0 0 1.0 0] zeros(2, 2)]

#animate_trajs_2d([[Ca * x for x in traj]], "LQ Tracking Trajectory", ["Agent", "Reference"], (-6, 6), (-6, 6), fps=20)

# Build a coordination sheaf for a multi-agent system.

# Build a sheaf with a circle topology
n_agents = 20
s = EuclideanSheaf{Float64}(repeat([4], n_agents))

for i in 2:n_agents
    add_sheaf_edge!(s, i - 1, i, C, C)
end
add_sheaf_edge!(s, 1, n_agents, C, C)

L = sheaf_laplacian_matrix(s)
global_state = BlockArray{Float64}(rand(-1.0:0.01:2.0, 4 * n_agents), repeat([4], n_agents))

function apply_L(L, x::BlockVector{Float64}, n_agents)
    return BlockArray{Float64}(L * x, repeat([4], n_agents))
end

function apply_L_nonlinear(L, x::BlockVector{Float64}, n_agents)
    return BlockArray{Float64}(L(x), repeat([4], n_agents))
end

# Set up dynamical system for iterating the sheaf Laplacian
γ = 1 / opnorm(L)
γ = γ - 0.1 * γ

function iterate_laplacian(L, x, n_agents, γ, niters)
    loss = [x' * L * x]
    for i in 1:niters
        x = x - γ * L * x
        push!(loss, x' * L * x)
    end
    return x, loss
end


function iterate_laplacian_nonlinear(L, x, n_agents, γ, niters)
    loss = [x' * L(x)]
    for i in 1:niters
        x = x - γ * L(x)
        push!(loss, x' * L(x))
    end
    return x, loss
end

#x_res, loss = iterate_laplacian(L, global_state, n_agents, γ, 200)

#plot(loss, yscale=:log10)

function run_sim1(L, x0, n_agents, niters, A, B, F)
    γ = 1 / opnorm(L)
    γ = γ - 0.1 * γ
    #x_target, _ = iterate_laplacian(L, x0, n_agents, γ, 200)
    traj = [x0]
    for i in 1:niters
        x_curr = traj[end]
        x_next = deepcopy(x_curr)
        x_target = x_curr - γ * apply_L(L, x_curr, n_agents)
        for i in 1:n_agents
            x_ref = x_target[Block(i)]
            x_ref[2] = 0.0
            x_ref[4] = 0.0
            u = F * (x_curr[Block(i)] - x_ref)
            x_next_i = A * x_curr[Block(i)] + B * u
            x_next[Block(i)] = x_next_i
        end
        push!(traj, x_next)
    end
    return traj
end

# Simulate communication delays.
function run_sim2(L, x0, n_agents, niters, A, B, F; comm_delay=10)
    γ = 1 / opnorm(L)
    γ = γ - 0.1 * γ
    #x_target, _ = iterate_laplacian(L, x0, n_agents, γ, 200)
    traj = [x0]
    x_target = x0 - γ * apply_L(L, x0, n_agents)
    for i in 1:niters
        x_curr = traj[end]
        x_next = deepcopy(x_curr)
        if i % comm_delay == 0
            x_target = x_curr - γ * apply_L(L, x_curr, n_agents)
        end
        for i in 1:n_agents
            x_ref = x_target[Block(i)]
            #=if i == 1
                x_ref[2] = 3.0
                x_ref[4] = 3.0
            else
                x_ref[2] = 0.0
                x_ref[4] = 0.0
            end=#

            x_ref[2] = 0.0
            x_ref[4] = 0.0
            u = F * (x_curr[Block(i)] - x_ref)
            x_next_i = A * x_curr[Block(i)] + B * u
            x_next[Block(i)] = x_next_i
        end
        push!(traj, x_next)
    end
    return traj
end

X, _ = ared(A, B, R, Q)
F = feedback_gain(A, B, R, X)
traj = run_sim2(L, global_state, n_agents, 1000, A, B, F)

t1 = [traj[i][Block(1)] for i in eachindex(traj)]
plt = plot_traj_2d(t1, C, "Agent Trajectories")
for i in 2:n_agents
    plot_traj_2d!(plt, [traj[j][Block(i)] for j in eachindex(traj)], C)
end
plt



# Flocking example with communication delays
function run_sim3(L, x0, n_agents, niters, A, B, F; comm_delay=10)
    #γ = 1 / opnorm(L)
    #γ = γ - 0.1 * γ
    γ = 0.0001
    x_target, loss = iterate_laplacian_nonlinear(L, x0, n_agents, γ, 10000)
    show(loss[end])
    traj = [x0]
    #x_target = x0 - γ * apply_L_nonlinear(L, x0, n_agents)
    for i in 1:niters
        x_curr = traj[end]
        x_next = deepcopy(x_curr)
        if i % comm_delay == 0
            x_target = x_target - γ * apply_L_nonlinear(L, x_curr, n_agents)
        end
        for i in 1:n_agents
            x_ref = x_target[Block(i)]
            #=if i == 1
                x_ref[2] = 1.0
                x_ref[4] = 1.0
            end=##=else
                                                                                                                            x_ref[2] = 0.0
                                                                                                                            x_ref[4] = 0.0
                                                                                                                        end=#

            #x_ref[2] = 0.0
            #x_ref[4] = 0.0
            u = F * (x_curr[Block(i)] - x_ref)
            x_next_i = A * x_curr[Block(i)] + B * u
            x_next[Block(i)] = x_next_i
        end
        push!(traj, x_next)
    end
    return traj
end

n_agents = 5
global_state = BlockArray{Float64}(rand(-2.0:0.01:2.0, 4 * n_agents), repeat([4], n_agents))

s = PotentialSheaf{EuclideanSheaf{Float64}}(repeat([4], n_agents))


q(x) = (x' * x - 5.0)^2
p(x) = [x[2], x[4]]' * [x[2], x[4]] + q([x[1], x[3]])
id = Matrix{Float64}(I(4))
for i in 2:n_agents
    add_sheaf_edge!(s, i - 1, i, id, id, p)
end
add_sheaf_edge!(s, 1, n_agents, id, id, p)

L = sheaf_laplacian(s)

traj = run_sim3(L, global_state, n_agents, 1000, A, B, F, comm_delay=1)


t1 = [traj[i][Block(1)] for i in eachindex(traj)]
plt = plot_traj_2d(t1, C, "Agent Trajectories")
for i in 2:n_agents
    plot_traj_2d!(plt, [traj[j][Block(i)] for j in eachindex(traj)], C)
end
plt

# Make a complete sheaf with barrier functions saying don't get too close
r = 5
U(y) = norm(y) < r ? exp(- 1 / (r^2 - norm(y)^2)) : 0.0
p(y) = 10*[y[2], y[4]]' * [y[2], y[4]] + U([y[1], y[3]])
V = [0 1.0 0 0; 0 0 0 1.0]

s = PotentialSheaf{EuclideanSheaf{Float64}}(repeat([4], n_agents))
for i in 1:n_agents
    for j in i+1:n_agents
        #=if i == 1
            add_sheaf_edge!(s, i, j, V, V, x-> 10* x'*x)
        else
            add_sheaf_edge!(s, i, j, C, C, U)
        end=#
        add_sheaf_edge!(s, i, j, id, id, p)
    end
end



#trajs = [[traj[i][Block(j)] for i in eachindex(traj)] for j in 1:n_agents]
#animate_trajs_2d(trajs, C, "Agent Trajectories", repeat([""], n_agents), (-2, 110), (-2, 110); fps=30)
#=
function compute_trajectory(A, B, F, x0, n; x_ref=zeros(length(x0)))
    traj = [x0]
    us = Vector{Float64}[]
    for i in 1:n-1
        x_curr = traj[end]
        u = F * (x_curr - x_ref)
        x_next = A * x_curr + B * u
        push!(traj, x_next)
        push!(us, u)
    end
    return traj, us
end

#traj, us = compute_trajectory(A, B, F, rand(4), 200, x_ref=[1.0, 0.0, 1.0, 0.0])
#plot_traj_2d(traj, C, "Single Agent LQR to (1,1)")


γ = 1 / opnorm(L)
γ = γ - 0.1 * γ
x_target, _ = iterate_laplacian(L, global_state, n_agents, γ, 200)

trajs = []

for i in 1:n_agents
    local x0 = global_state[Block(i)]
    x_ref = x_target[Block(i)]
    x_ref[2] = 0.0
    x_ref[4] = 0.0
    local traj, _ = compute_trajectory(A, B, F, x0, 200, x_ref=x_ref)
    push!(trajs, traj)
end

plt = plot_traj_2d(trajs[1], C, "Multi-Agent LQR to Sheaf Consensus")
for i in 2:n_agents
    plot_traj_2d!(plt, trajs[i], C)
end
plt=#

# Try a proximal method

function augmented_A(A, N)
    res = hcat(I(4), zeros(4, 4(N-1)))
    for i in 2:N
        row = hcat(zeros(4, 4(i-2)), -A, I(4), zeros(4, 4(N-i)))
        res = vcat(res, row)
    end
    return res
end

function augmented_B(B, N)
    return Array(sparse(blocksparse(collect(1:N), collect(1:N), repeat([-B], N))))
end

function dynamic_constraint(A, B, N, x0)
    return [augmented_A(A, N) augmented_B(B, N)], vcat(A*x0, zeros(4 * (N-1)))
end
#c = rand(4)
#C, b = dynamic_constraint(A, B, 5, c)

#ind = IndAffine(C, b)
#x0 = rand(30)
#r, _ = prox(ind, x0)

function run_sim4(A,B,N,s,x0,niters,n_agents)
    L = sheaf_laplacian(s)
    γ = 0.01
    state_blocks = repeat([4], n_agents)

    traj = [x0]
    for i in 1:niters
        x_curr = deepcopy(traj[end])
        x_target = x_curr - γ * apply_L_nonlinear(L, x_curr, n_agents)

        # Make these blocked for the agents to access
        x_curr_b = BlockArray{Float64}(x_curr, state_blocks)
        x_target_b = BlockArray{Float64}(x_target, state_blocks)

        for j in 1:n_agents
            x_ref = x_target_b[Block(j)]
            if j == 1
                x_ref[2] = 1.0
                x_ref[4] = 0.0
            end
            x_state = x_curr_b[Block(j)]

            # Build the dynamics constraint using x_state as x0
            C, b = dynamic_constraint(A, B, N, x_state)
            ind = IndAffine(C, b)

            # Prepare input for prox
            target = vcat(repeat([x_ref], N)..., zeros(2* N))
            r, _ = prox(ind, target)

            # Get the first control input
            u = r[4N+1:4N+2]
            x_next = A * x_state + B * u
            x_curr_b[Block(j)] = x_next
        end
        push!(traj, x_curr_b)

    end
    return traj
end



traj = run_sim4(A,B,10,s,global_state,1000,n_agents)


t1 = [traj[i][Block(1)] for i in eachindex(traj)]
plt = plot_traj_2d(t1, C, "Agent Trajectories")
for i in 2:n_agents
    plot_traj_2d!(plt, [traj[j][Block(i)] for j in eachindex(traj)], C)
end
plt


trajs = [[traj[i][Block(j)] for i in eachindex(traj)] for j in 1:n_agents]
animate_trajs_2d(trajs, C, "Agent Trajectories", repeat([""], n_agents), (-6, 100), (-6, 6); fps=30)

#=
function augmented_A(A, N)
    return vcat(I(4), [A^n for n in 1:N]...)
end

function augmented_B(A, B, N)
    column(i) = vcat(zeros(4 * (i), size(B, 2)), [A^(j - i) * B for j in i:N]...)
    return hcat([column(i) for i in 1:N]...)
end=#