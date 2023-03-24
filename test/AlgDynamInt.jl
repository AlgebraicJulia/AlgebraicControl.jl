using Zygote
using AlgebraicDynamics
using AlgebraicDynamics.UWDDynam
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Programs


A = [0 1; .01 0]
B = I(2)
Q = [1.0 0; 0 1.0]
R = [1.0]
x₀ = [3.0, 0]

dotx(x,u,t) = A*x + B*u[1:2]
doty(y,u,t) = A*y + B*u[3:4]

s1 = ContinuousResourceSharer{Float64}(2, dotx)
s2 = ContinuousResourceSharer{Float64}(2, doty)


composition_pattern = @relation (x, z) begin
    s1(x,y)
    s2(y,z)
end

submodels = Dict(:s1 => s1, :s2 => s2)

to_graphviz(composition_pattern, box_labels = :name, junction_labels = :variable, edge_attrs=Dict(:len => ".75"))

system = oapply(composition_pattern, submodels)

systemd = euler_approx(system, 0.1)

x1 = eval_dynamics(systemd, ones(3), ones(4))

#g = gradient(x -> eval_dynamics(systemd,x,ones(4)), ones(3))

s = Box(:sensor, [:s, :e], [:s′])
c = Box(:controller, [:d, :s′], [:c])
d = Box(:dynamics, [:c], [:s]);

# A wiring diagram has outer inports and outports which define the interface of target system. 
# Then we add the boxes and wires to the diagram and visualize the result.

UAV = WiringDiagram([:e,:d], [:s])

sensor     = add_box!(UAV, s)
controller = add_box!(UAV, c)
dynamics   = add_box!(UAV, d)

add_wires!(UAV, [
    ## net inputs
    (input_id(UAV), 1) => (sensor, 2),
    (input_id(UAV), 2) => (controller, 2),

    ## connections
    (sensor, 1) => (controller, 1),
    (controller, 1) => (dynamics, 1),
    (dynamics, 1) => (sensor, 1),

    ## net output
    (dynamics, 1) => (output_id(UAV), 1)
]);

#-

to_graphviz(UAV)


function 𝗟(𝐖)
    𝐿(u, x, p, t) = [ -p.𝓐l * (u[1] - x[1] - x[2]) ] # sc
    𝐶(u, x, p, t) = [ -p.𝓐c * (u[1] + p.𝓑c*x[1] - x[2]) ] # sl
    𝐷(u, x, p, t) = LVector(α = -0.313*u[1] +  56.7*u[2] +  0.232*x[1],
                             q = -0.013*u[1] - 0.426*u[2] + 0.0203*x[1],
                             θ =  56.7*u[2]              )

    u_𝐿(u,p,t) = [ u[1] ] # outputs sl
    u_𝐶(u,p,t) = [ u[1] ] # outputs sc
    u_𝐷(u,p,t) = [ u[3] ] # outputs θ

    return oapply(𝐖,
                  Dict(:sensor     => ContinuousMachine{Float64}(2, 1, 1, 𝐿, u_𝐿),
                       :controller => ContinuousMachine{Float64}(2, 1, 1, 𝐶, u_𝐶),
                       :dynamics   => ContinuousMachine{Float64}(1, 3, 1, 𝐷, u_𝐷)))
end

uav_system = 𝗟(UAV)


