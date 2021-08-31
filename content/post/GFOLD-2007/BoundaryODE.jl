using DifferentialEquations, Plots, LinearAlgebra

tspan = (0.0,1.0) # converted to fixed time problem by treating final time as ...

g = [0, 0, -3.7114];
cosphi = cos(deg2rad(27));
α = 1 /(225 * 9.807 * cosphi);

T_max = 3100;
n = 6;
T1 = 0.3 * T_max;
T2 = 0.8 * T_max;
ρ₁ = n * T1 * cosphi;
ρ₂ = n * T2 * cosphi;

mwet = 1905;

r₀ = [2000, 0, 1500] # x y z, unlike in paper BE CAREFUL
rdot₀ = [100, 0, -75]

rfinal = [0, 0, 0]
rdotfinal = [0, 0, 0]

λ₀ = -1
lambda_guess = randn(7) 

t₀ = 0
Δt_guess = 72 # Nonlinear Optimal Control Theory by Leonard David Berkovitz, Negash G. Medhin (z-lib.org).pdf pg 24

function simplependulum!(du,u,p,s)
    r = u[1:3]
    rdot = u[4:6]
    m = u[7]
    λ₁ = u[8:10]
    λ₂ = u[11:13]
    λ₃ = u[14]
    Δt = u[15]
    real_t = t₀ + Δt * s

    (λ₀, α, ρ₁, ρ₂) = p

    R12 = λ₀ - α * λ₃ + norm(λ₂ / m);
    Γ = if (R12 > 0) ρ₂ else ρ₁ end
    T_c = λ₂ / norm(λ₂) * Γ

    du[1:3] = rdot
    du[4:6] = g + T_c / m
    du[7] = -α * Γ
    du[8:10] = [0 0 0]
    du[11:13] = -λ₁
    du[14] = dot(λ₂, T_c) / m^2 
    du[15] = 0.0

    du *= Δt
end

function bc2!(residual, u, p, t) # u[1] is the beginning of the time span, and u[end] is the ending
    residual[1:3] = u[1][1:3] - r₀
    residual[4:6] = u[1][4:6] - rdot₀
    residual[7] = u[1][7] - mwet
    residual[8:10] = u[end][1:3] - rfinal
    residual[11:13] = u[end][4:6] - rdotfinal
    residual[14] = u[end][14] - 0
    residual[15] = if (u[1][15] > 0) 0 else u[1][15] end
    #residual[16] = u[end][15] - u[1][15] # IDK WHY THIS ISN'T WORKING
end

bvp = BVProblem(simplependulum!, bc2!, [r₀; rdot₀; mwet; lambda_guess; Δt_guess], tspan, (λ₀, α, ρ₁, ρ₂))
sol = solve(bvp, GeneralMIRK4(), dt=0.01)

plotlyjs()
plot(sol)