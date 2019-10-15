using PendulumDyn

# Simulation parameters
dt = 0.01 # step size
t = 0:dt:10 # simulation steps
u(t) = 10.0*sin(5*t) # input (voltage to motor)

# Create a pendulum object, contains standard parameters for rod length, mass etc
pend = Pendulum(Float64)

# Preallocate memory to store states
X = Array{Float64}(undef, 4, length(t))
xt = zeros(Float64, 4)

# Simulate
for i = 1:length(t)
    X[:, i] = xt
    pendulum_rk4!(xt, u(t[i]), dt, pend)
end

# Plot results
using Plots

fig = plot(layout = (2, 1))
plot!(fig[1], t, X[1, :], lab = "position")
plot!(fig[2], t, X[3, :], lab = "angle")
