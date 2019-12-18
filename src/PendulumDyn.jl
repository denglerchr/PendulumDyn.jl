module PendulumDyn

include("pendulummotor.jl")
export dxdt_pendulummotor, pendulummotor_rk4, pendulummotor_rk4!, PendulumMotor

include("pendulum.jl")
export dxdt_pendulum, dxdt_pendulum!, pendulum_rk4, pendulum_rk4!, Pendulum

include("ode_sol.jl")
export pendulum_odesol, pendulum_odesol!

end
