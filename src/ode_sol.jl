# Do multiple small timesteps in case of a large dt
function pendulum_odesol!(xt::AbstractVector, ut::Number, dt::Number, pend::Pendulum{T} = Pendulum()) where T
    if dt < T(0.06) # only one rk4 step
        pendulum_rk4!(xt, ut, dt, pend)
    else # multiple one rk4 step
        N = div(dt, T(0.06)) + one(T)
        dt2 = dt/N
        for j = 1:N
            pendulum_rk4!(xt, ut, dt2, pend)
        end
    end
    return xt
end

function pendulum_odesol(xt::AbstractVector, ut::Number, dt::Number, pend::Pendulum{T} = Pendulum()) where T
    xt2 = copy(xt)
    return pendulum_odesol!(xt2, ut, dt, pend)
end

#=
using OrdinaryDiffEq


# Use a different ODE solver
function pendulum_odesol2(xt::AbstractVector{T}, ut::T, dt::T, pend::Pendulum{T} = Pendulum(); alg = DP5()) where {T<:Real}
    # Define ode problem
    tspan = (zero(T), dt)
    dynfunc(x, p, t) = dxdt_pendulum(x, ut, pend)
    prob = ODEProblem(dynfunc, xt, tspan)
    sol = solve(prob, alg ,reltol=1e-8,abstol=1e-8)

    #return last state
    return sol.u[end]
end
=#
