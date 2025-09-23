using AlgebraicOptimization
using AlgebraicControl
using Convex
using LinearAlgebra
using ProximalAlgorithms
using MatrixEquations
using SparseArrays
using Plots
using BlockArrays

# Set up each agent's dynamics: x' = Ax + Bu
dt = 0.1  # Discretization step size
A_sub = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]
B_sub = [0 0; dt 0; 0 0; 0 dt]
C = [1.0 0 0 0; 0 0 1.0 0] # Output agent's position

function block_diag(A, n)
    return Array(sparse(blocksparse(collect(1:n), collect(1:n), repeat([A], n))))
end

A = block_diag(A_sub, 3)
B = block_diag(B_sub, 3)


s = @cellular_sheaf C begin
    x::Stalk{4}, y::Stalk{4}, z::Stalk{4}

    C(x) == C(y)
    C(x) == C(z)
    C(y) == C(z)
end

L = Array(sparse(sheaf_laplacian_matrix(s)))
V_sub = I(4)
V_sub[1, 1] = 0.0
V_sub[3, 3] = 0.0
V = block_diag(V_sub, 3)

Q = L + V

R = I(6)

X, _ = ared(A, B, R, Q)

F = -inv(R + B' * X * B) * (B' * X * A)

function sparsity_pattern(A_in, zero_threshold)
    A = deepcopy(A_in)
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            if abs(A[i, j] < zero_threshold)
                A[i, j] = 0.0
            end
        end
    end
    return sparse(A)
end

function compute_trajectory(A, B, F, x0, n)
    traj = [x0]
    us = Vector{Float64}[]
    for i in 1:n-1
        x_curr = traj[end]
        u = F * x_curr
        x_next = A * x_curr + B * u
        push!(traj, x_next)
        push!(us, u)
    end
    return traj, us
end

traj, us = compute_trajectory(A, B, F, rand(-2.0:0.1:2.0, 12), 100)

# Postprocess for plotting
traj1 = [C * x[1:4] for x in traj]
traj2 = [C * x[5:8] for x in traj]
traj3 = [C * x[9:12] for x in traj]

traj1 = mapreduce(permutedims, vcat, traj1)
traj2 = mapreduce(permutedims, vcat, traj2)
traj3 = mapreduce(permutedims, vcat, traj3)

plt1 = plot(traj1[:, 1], traj1[:, 2])
plot!(plt1, traj2[:, 1], traj2[:, 2])
plot!(plt1, traj3[:, 1], traj3[:, 2])

# Build a sheaf with a circle topology
n_agents = 21
s = EuclideanSheaf{Float64}(repeat([4], n_agents))

for i in 2:n_agents
    add_sheaf_edge!(s, i - 1, i, C, C)
end
add_sheaf_edge!(s, 1, n_agents, C, C)

L = Array(sparse(sheaf_laplacian_matrix(s)))
V = block_diag(V_sub, n_agents)
Q = L + V
R = I(n_agents * 2)

A = block_diag(A_sub, n_agents)
B = block_diag(B_sub, n_agents)

X, _ = ared(A, B, R, Q)

F = -inv(R + B' * X * B) * (B' * X * A)

#F = sparsity_pattern(F, 1e-8)


Fb = BlockMatrix(F, repeat([2], n_agents), repeat([4], n_agents))
Fb_adj = deepcopy(Fb)
for i in 1:n_agents
    for j in 1:n_agents
        if i == j || (i + 1) % n_agents == j || (i - 1) % n_agents == j #|| (i + 2) % n_agents == j || (i - 2) % n_agents == j
            continue
        end
        Fb_adj[Block(i), Block(j)] = zeros(2, 4)
    end
end

Fb_adj_twohop = deepcopy(Fb)
for i in 1:n_agents
    for j in 1:n_agents
        if i == j || (i + 1) % n_agents == j || (i - 1) % n_agents == j || (i + 2) % n_agents == j || (i - 2) % n_agents == j
            continue
        end
        Fb_adj_twohop[Block(i), Block(j)] = zeros(2, 4)
    end
end

Fb_diag = deepcopy(Fb)
for i in 1:n_agents
    for j in 1:n_agents
        if i == j
            continue
        end
        Fb_diag[Block(i), Block(j)] = zeros(2, 4)
    end
end

x0 = rand(-2.0:0.1:4.0, n_agents * 4)
traj_decentralized, us_decentralized = compute_trajectory(A, B, Fb_diag, x0, 400)
traj_onehop, us_onehop = compute_trajectory(A, B, Fb_adj, x0, 400)
traj_twohop, us_twohop = compute_trajectory(A, B, Fb_adj_twohop, x0, 400)
function plot_traj(traj, n_agents, title)
    traj_mat = mapreduce(permutedims, vcat, traj)
    slices = []
    for i in 1:n_agents
        push!(slices, (i-1)*4+1:i*4)
    end
    agent_trajs = [mapreduce(permutedims, vcat, [C * x for x in eachrow(traj_mat[:, slice])]) for slice in slices]

    plt2 = plot(title=title)
    for traj in agent_trajs
        plot!(plt2, traj[:, 1], traj[:, 2], label="")
    end

    return plt2
end

function coordination_loss(Q, traj)
    return sum([x[1]' * Q * x[1] for x in eachrow(traj)])
end

function controls_loss(R, us)
    return sum([u' * R * u for u in us])
end

plt_decentralized = plot_traj(traj_decentralized, n_agents, "Decentralized")
plt_onehop = plot_traj(traj_onehop, n_agents, "1-Hop")
plt_twohop = plot_traj(traj_twohop, n_agents, "2-Hop")

traj_centralized, us_centralized = compute_trajectory(A, B, F, x0, 400)

centralized_loss = coordination_loss(Q, traj_centralized) + controls_loss(R, us_centralized)
decentralized_loss = coordination_loss(Q, traj_decentralized) + controls_loss(R, us_decentralized)

println("Centralized loss: $centralized_loss")
println("Decentralized loss: $decentralized_loss")

plt_centralized = plot_traj(traj_centralized, n_agents, "Centralized")

#plot(plt_centralized, plt_decentralized, plt_onehop, plt_twohop, layout=@layout [a b; c d])

#=Xb = BlockMatrix(X, repeat([4], n_agents), repeat([4], n_agents))
Xb_adj = deepcopy(Xb)

for i in 1:n_agents
    for j in 1:n_agents
        if i == j || (i + 1) % n_agents == j || (i - 1) % n_agents == j || (j + 1) % n_agents == i || (j - 1) % n_agents == i
            continue
        end
        Xb_adj[Block(i), Block(j)] = zeros(4, 4)
    end
end=#

function is_dd(G::Matrix{Float64})
    return all(sum(view(G, i, :)) <= 2abs(G[i, i]) for i in axes(G, 1))
end

function rotation_matrix(theta)
    return [cos(theta) -sin(theta); sin(theta) cos(theta)]
end

n_agents = 5
# Build a sheaf with rotation matrices on edges
function make_rotation_sheaf(n_agents)
    s = EuclideanSheaf{Float64}(repeat([4], n_agents))
    θ = 2π / n_agents
    θ0 = 0.0

    for i in 2:n_agents
        θ0 += θ
        r = rotation_matrix(θ0)
        add_sheaf_edge!(s, i - 1, i, C, r * C)
    end
    add_sheaf_edge!(s, 1, n_agents, C, rotation_matrix(θ0 - θ) * C)
    return s
end

s = make_rotation_sheaf(n_agents)

L = Array(sparse(sheaf_laplacian_matrix(s)))
V = block_diag(V_sub, n_agents)
Q = L #+ V
R = I(n_agents * 2)

A = block_diag(A_sub, n_agents)
B = block_diag(B_sub, n_agents)

X, _ = ared(A, B, R, Q)

F = -inv(R + B' * X * B) * (B' * X * A)

traj, us = compute_trajectory(A, B, F, rand(-2.0:0.1:2.0, n_agents * 4), 100)

plot_traj(traj, n_agents, "Rotation Sheaf")

# Example from Matrix Weighted Consensus paper

A12 = [2 0; 0 1]
A13 = [2 3; 3 5]
A47 = [0 0; 0 1]
A56 = [1 0; 0 0]
A14 = [0.75 -0.433; -0.433 0.25]
A17 = [0.75 0.433; 0.433 0.25]
A45 = [1 0.5; 0.5 1]
A46 = [0.9518 -0.2142; -0.2142 0.0482]
A78 = [3 2; 2 3]
A89 = [2 0; 0 2]

s = EuclideanSheaf{Float64}(repeat([4], 9))

function add_matrix_weighted_edge!(s, i, j, A)
    _, U = qr(A)
    add_sheaf_edge!(s, i, j, U * C, U * C)
end

add_matrix_weighted_edge!(s, 1, 2, A12)
add_matrix_weighted_edge!(s, 1, 3, A13)
add_matrix_weighted_edge!(s, 4, 7, A47)
add_matrix_weighted_edge!(s, 5, 6, A56)
add_matrix_weighted_edge!(s, 1, 4, A14)
add_matrix_weighted_edge!(s, 1, 7, A17)
add_matrix_weighted_edge!(s, 4, 5, A45)
add_matrix_weighted_edge!(s, 4, 6, A46)
add_matrix_weighted_edge!(s, 7, 8, A78)
add_matrix_weighted_edge!(s, 8, 9, A89)

n_agents = 9
L = Array(sparse(sheaf_laplacian_matrix(s)))
V = block_diag(V_sub, n_agents)
Q = L + V
R = .3 .* I(n_agents * 2)

A = block_diag(A_sub, n_agents)
B = block_diag(B_sub, n_agents)

X, _ = ared(A, B, R, Q)

F = -inv(R + B' * X * B) * (B' * X * A)

traj, us = compute_trajectory(A, B, F, rand(-2.0:0.1:2.0, n_agents * 4), 400)

plt = plot_traj(traj, n_agents, "Matrix Weighted Consensus Sheaf")


function compute_axis_limits(agent_trajs; margin=1.0)
    xs = vcat([traj[:, 1] for traj in agent_trajs]...)
    ys = vcat([traj[:, 3] for traj in agent_trajs]...)
    x_min, x_max = minimum(xs), maximum(xs)
    y_min, y_max = minimum(ys), maximum(ys)
    return (x_min - margin, x_max + margin), (y_min - margin, y_max + margin)
end

function animate_trajectory(trajectory, n_agents; fps=20)
    colors = [:red, :blue, :green, :orange, :purple, :black, :magenta, :cyan, :brown, :gray]
    traj_mat = mapreduce(permutedims, vcat, trajectory)
    agent_trajs = [traj_mat[:, (i-1)*4 .+ (1:4)] for i in 1:n_agents]
    xlims, ylims = compute_axis_limits(agent_trajs)
    anim = @animate for t in 1:length(trajectory)
        plt = plot(
            title="Agent Consensus Over Time",
            xlabel="x",
            ylabel="y",
            legend=false,
            xlims=xlims,
            ylims=ylims,
        )
        for i in 1:n_agents
            pos = agent_trajs[i][t, :]
            scatter!(
                plt,
                [pos[1]], [pos[3]],
                color=colors[mod1(i, length(colors))],
                ms=5,
            )
        end
    end

    return gif(anim,"anim.gif", fps=fps)
end

#animate_trajectory(traj, n_agents, fps=20)

animate_trajectory(traj_centralized[1:100], 21, fps=20)


