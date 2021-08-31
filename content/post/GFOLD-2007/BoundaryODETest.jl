using DifferentialEquations, Plots, LinearAlgebra

tspan = (0.0,1.0) # converted to fixed time problem by treating final time as ...

t₀ = 0
Δt_guess = 5 # Nonlinear Optimal Control Theory by Leonard David Berkovitz, Negash G. Medhin (z-lib.org).pdf pg 24

function simplependulum!(du,u,p,s)
    x₁, x₂, λ₁, λ₂, Δt = u
    real_t = t₀ + Δt * s

    du = [x₂, -λ₂, 0, -λ₁, 0]
    du *= Δt
end

function bc2!(residual, u, p, t) # u[1] is the beginning of the time span, and u[end] is the ending
    residual[1] = u[1][1] - 1.0
    residual[2] = u[1][2] - 2.0
    residual[3] = u[end][1] - 3.0
    residual[4] = u[end][4] - 0.0
end

bvp = BVProblem(simplependulum!, bc2!, [1; 2; randn(2); Δt_guess], tspan)
sol = solve(bvp, GeneralMIRK4(), dt = 0.01)

plotlyjs()
plot(sol)