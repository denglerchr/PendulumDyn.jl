struct PendulumMotor{T}
    # Constant mechanical parameters
    Mcart::T # mass of the cart
    Mrod::T # mass of the rod
    Lrod::T    # Length of the rod
    mu_v::T  # viscous friction

    # Constant electric parameters
    k1::T # motor constant
    k2::T # motor constant
    i::T # gear ratio
    r::T
    R::T # armature resistance
    L::T # armature inductance
end

function PendulumMotor(T::Type;
    Mcart = 29.691,  # mass of the cart
    Mrod = 2.3*1.178,  # mass of the rod
    Lrod = 1.178/2,   # Length of the rod
    mu_v = 0.001, # viscous friction
    k1 = 0.263, # motor constant
    k2 = 0.268, # motor constant
    i = 3, # gear ratio
    r = 0.0239,
    R = 1.35, # armature resistance
    L = 1.35*10^-3 ) # armature inductance
    return PendulumMotor{T}(T(Mcart), T(Mrod), T(Lrod), T(mu_v), T(k1), T(k2), T(i), T(r), T(R), T(L))
end

PendulumMotor() = PendulumMotor(Float64)


function dxdt_pendulummotor(state::AbstractVector{T}, u::T) where {T}
    pend = PendulumMotor(T)
    return dxdt_pendulum(state, u, pend)
end

function dxdt_pendulummotor!(state::AbstractVector{T}, u::T) where {T}
    pend = PendulumMotor(T)
    return dxdt_pendulum!(state, u, pend)
end

function dxdt_pendulum(state::AbstractVector, u::T, pend::PendulumMotor{T}) where {T}
    state2 = copy(state)
    return  dxdt_pendulum!(state2, u, pend)
end

function dxdt_pendulum!(state::AbstractVector, u::T, pend::PendulumMotor{T}) where {T}
    g = T(9.81)

    # Rename some derivatives
    x = state[1]
    dx = state[2]
    phi = state[3]
    dphi = state[4]
    current = state[5]
    u2 = clamp.(u, -T(40), T(40))

    # Motor to pend
    Mm = pend.k2*current # k2*i is motor moment
    Fout = pend.i/pend.r*Mm # force by motor on the pend

    # Pendulum dynamics (from "http://sun4.vaniercollege.qc.ca/~iti/proj/Denis-cart-pendulum-spring.pdf")
    K = (pend.Mcart + pend.Mrod*sin(phi)^2)
    dv = (pend.Mrod*pend.Lrod*dphi^2*sin(phi) + pend.Mrod*g*sin(phi)*cos(phi) - pend.mu_v*dx + Fout)/K
    ddphi = (-pend.Mrod*pend.Lrod*dphi^2*sin(phi)*cos(phi) - (pend.Mcart + pend.Mrod)*g*sin(phi) + (pend.mu_v*dx - Fout)*cos(phi) )/(K*pend.Lrod)

    # Current dynamics
    omega = dx*pend.i/pend.r
    di::T = -pend.R/pend.L * current - pend.k1/pend.L * omega + u2/pend.L

    state[1] = dx
    state[2] = dv
    state[3] = dphi
    state[4] = ddphi
    state[5] = di
    return state
end

function pendulummotor_rk4(xt::AbstractVector, ut::T, dt::T) where T
    pend = PendulumMotor(T)
    return pendulum_rk4(xt, ut, dt, pend)
end

function pendulummotor_rk4!(xt::AbstractVector, ut::T, dt::T) where T
    pend = PendulumMotor(T)
    return pendulum_rk4!(xt, ut, dt, pend)
end

function pendulum_rk4(xt::AbstractVector, ut::Number, dt::Number, pend::PendulumMotor{T}) where T
    xt2 = copy(xt)
    return pendulum_rk4!(xt2, ut, dt, pend)
end

function pendulum_rk4!(xt::AbstractVector, ut::Number, dt::Number, pend::PendulumMotor{T}) where T
    temp = similar(xt)

    #get partial results dy1-dy4
    dy1 = dxdt_pendulum(xt, ut, pend);
    temp .= dy1.*(dt/2) .+ xt;

    dy2 = dxdt_pendulum(temp, ut, pend);
    temp .= dy2.*(dt/2) .+ xt;

    dy3 = dxdt_pendulum(temp, ut, pend);
    temp .= dy3.*dt .+ xt;

    dy4 = dxdt_pendulum(temp, ut, pend);

    # update xt
    @. xt += dt*( dy1 + 2.0*(dy2+dy3)+dy4 )/6.0
    return xt
end
