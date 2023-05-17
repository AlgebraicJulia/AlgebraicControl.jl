using Plots
import Plots: plot
using LinearAlgebra
import Base: *, +

struct MPCBifunction
    state_space::Int
    control_space::Int
    cost::Function # state_space × control_space → R
    dynamics::Function # state_space × control_space → state_space
end

struct Pack
    nports::Int # P
    lb::Vector{Float64} # lb : P → R
    ub::Vector{Float64} # ub : P → R
    step_sizes::Vector{Float64}
    resolution::Vector{Int} # resolution : P → N
end
function Pack(nports::Int, lb::Vector{Float64}, ub::Vector{Float64}, resolution::Vector{Int})
    step_sizes = Float64[]
    for i in 1:nports
        push!(step_sizes, (ub[i] - lb[i]) / resolution[i])
    end
    return Pack(nports, lb, ub, step_sizes, resolution)
end

function entries(p::Pack)
    ranges = [1:u for u in p.resolution]
    return Iterators.product(ranges...)
end


struct Pixel
    bounds::Vector{Pair{Float64}}
end

# Get the pixel of an entry as a list of lb=>ub pairs for each dimension
function Pixel(p::Pack, e::Vector{Int})
    result = Pair{Float64}[]
    for i in 1:p.nports
        push!(result, (p.lb[i] + p.step_sizes[i]*(e[i]-1)) => p.lb[i] + p.step_sizes[i]*e[i])
    end
    return Pixel(result)
end

+(p1::Pair{Float64}, p2::Pair{Float64}) = (p1[1] + p2[1]) => (p1[2] + p2[2])

function evaluate_pixel(p::Pixel, f::MPCBifunction; ϵ=0.1)
    # Get centroid of pixel
    centroid = Float64[]
    for (lb,ub) in p.bounds
        push!(centroid, (lb+ub)/2)
    end
    #println(centroid)
    # Group centroid into inputs for the bifunction
    x_in = centroid[1:f.state_space]
    u = centroid[f.state_space+1:f.state_space+f.control_space]
    x_out_desired = centroid[f.state_space+f.control_space+1:end]
    x_out = f.dynamics(x_in, u)
    if norm(x_out - x_out_desired) <= ϵ
        return f.cost(x_in, u)
    else
        return Inf
    end
end

struct PackUWD
    dom::Vector{Pack}
    codom::Pack
    links::Pack
    φ::Function
end

# Super naive implementations for now!
struct PixelMatrix
    x_res::Int
    y_res::Int
    x_range::Pair{Float64}
    y_range::Pair{Float64}
    m::Matrix{Bool}
end

function PixelMatrix(x_res, y_res, x_range, y_range, rel::Function)
    m = repeat(repeat([false], x_res)', y_res)
    x_step = (x_range[2] - x_range[1]) / x_res
    y_step = (y_range[2] - y_range[1]) / y_res
    for (x,i) in zip(range(x_range[1], x_range[2], x_res), 1:x_res)
        for (y,j) in zip(range(y_range[1], y_range[2], y_res), 1:y_res)
            s1 = rel(x,y) < 0
            s2 = rel(x+x_step,y) < 0
            s3 = rel(x,y+y_step) < 0
            s4 = rel(x+x_step,y+y_step) < 0
            #println("Checking the relation on ($x,$y).")
            #m[i,j] = rel(x,y)
            #m[i,j] = s1 != s4
            m[i,j] = s1 != s4 || s2 != s3
            #m[i,j] = s1 != s2 || s3 != s4
            #m[i,j] = s1 != s3 || s2 != s4
        end
    end
    return PixelMatrix(x_res, y_res, x_range, y_range, m)
end

function *(pm1::PixelMatrix, pm2::PixelMatrix)
    pm1.y_res == pm2.x_res || error("Dimension mismatch")
    m = pm1.m*pm2.m .≠ 0
    return PixelMatrix(pm1.x_res, pm2.y_res, pm1.x_range, pm2.y_range, m)
end

function plot(pm::PixelMatrix)
    xs = Int[]
    ys = Int[]

    for i in 1:size(pm.m)[1]
        for j in 1:size(pm.m)[2]
            if pm.m[i,j]
                push!(xs, i)
                push!(ys, j)
            end
        end
    end
    return scatter(xs,ys)
end

# Example usage
r1(x,w) = x^2 - w
r2(w,y) = 1 - y^2 - w

pm1 = PixelMatrix(100,100,-1.5=>1.5,-1.5=>1.5,r1)
pm2 = PixelMatrix(100,100,-1.5=>1.5,-1.5=>1.5,r2)

res = pm1*pm2
plot(res)

# Example bifunction usage
f = MPCBifunction(1,1,(x,u)->x[1]^2 + u[1]^2, (x,u)->[5*x[1] + u[1]])
p = Pack(3, repeat([-1.5],3), repeat([1.5],3), repeat([100],3))

for e in entries(p)
    pix = Pixel(p, collect(e))
    val = evaluate_pixel(pix, f, ϵ=0.015)
    if val != Inf
        println(pix)
        println(val)
    end
end