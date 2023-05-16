using Plots
import Plots: plot
import Base: *

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