using DifferentialEquations, Plots, LinearAlgebra, XLSX

tspan = (0.0, 3.45)

g = [0, 0, -9.80655];

ThrustAcceleration = 14

r₀ = [5, -4, 20] # x y z, unlike in paper BE CAREFUL
rdot₀ = [2, 1, -10]

rfinal = [0, 0, 0]
rdotfinal = [0, 0, 0]

lambda_guess = randn(6) 

function simplependulum!(du,u,p,s)
    r = u[1:3]
    rdot = u[4:6]
    λ₁ = u[7:9]
    λ₂ = u[10:12]

    control = λ₂ / norm(λ₂) * ThrustAcceleration

    du[1:3] = rdot
    du[4:6] = g + control
    du[7:9] = [0 0 0]
    du[10:12] = -λ₁
end

function bc2!(residual, u, p, t) # u[1] is the beginning of the time span, and u[end] is the ending
    residual[1:3] = u[1][1:3] - r₀
    residual[4:6] = u[1][4:6] - rdot₀
    #cost = -1 * norm(u[end][4:6])
    # cost = -1 * 2 * dot(u[1][4:6], (g + u[1][10:12] / norm(u[1][10:12]) * ThrustAcceleration))
    #residual[7] = cost + dot(u[end][7:9], u[end][4:6]) + dot(u[end][10:12], (g + u[end][10:12] / norm(u[end][10:12]) * ThrustAcceleration))
    residual[7:9] = u[end][4:6] - rdotfinal
    residual[10:12] = u[end][1:3] - rfinal
end

bvp = BVProblem(simplependulum!, bc2!, [r₀; rdot₀; lambda_guess], tspan)
sol = solve(bvp, GeneralMIRK4(), dt=0.4)

plotlyjs()
u = transpose(reduce(hcat, sol.u))
control = u[:, 10:12] ./ norm.(eachrow(u[:, 10:12])) #* ThrustAcceleration
plot(sol.t, norm.(eachrow(control)))
plot(sol.t, u[:, 4:6])