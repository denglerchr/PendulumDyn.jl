using PendulumDyn

# Simulation parameters
dt = 0.01 # step size
Ntimes = 1

dt2 = Ntimes*dt
t = 0:dt:10 # simulation steps
t2 = 0:dt2:10 # simulation steps

u(t) = 10.0*sin(5*t) # input (voltage to motor)

# Create a pendulum object, contains standard parameters for rod length, mass etc
pend = Pendulum(Float64)

# Preallocate memory to store states
X = Array{Float64}(undef, 4, length(t))
X2 = Array{Float64}(undef, 4, length(t2))

xt = zeros(Float64, 4)
xt2 = copy(xt)

# Simulate
@btime for i = 1:length(t)
    X[:, i] = xt
    pendulum_rk4!(xt, u(t[i]), dt, pend)
end

@btime for i = 1:length(t2)
    X2[:, i] = xt2
    pendulum_odesol!(xt2, u(t2[i]), dt2, pend)
end

# Plot results
using Plots

fig = plot(layout = (2, 1))
plot!(fig[1], t, X[1, :], lab = "position")
plot!(fig[1], t2, X2[1, :], lab = "")

plot!(fig[2], t, X[3, :], lab = "angle")
plot!(fig[2], t2, X2[3, :], lab = "")
