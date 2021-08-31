using DifferentialEquations

x0 = [1.0, 2.0, 0.0]

t_start = 0.0
t_final = 10.0
tspan = (t_start, t_final)

u = t -> sin(t)

function f(xdot, x, u, t)
    xdot[1] = 0.01 * x[1] - u(t)
    xdot[2] = -.05 * x[2] + u(t)
    xdot[3] = transpose(x[1:end-1]) * x[1:end-1] + transpose(u(t)) * u(t)
end

prob = ODEProblem(f, x0, tspan, u)

p = Iterators.product(-5.0:0.25:5.0, -5.0:0.25:5.0, 0.0);
x0_range = vec(collect.(p))

u0_range = range(-5.0, 5.0, length = 101)

numberOfParameters = length(x0_range) * length(u0_range)
print(numberOfParameters)

function parameterChange(prob, i, repeat)
    remake(prob, u0=x0_range[(i - 1) % length(x0_range) + 1], p = t -> u0_range[1 + (i - 1) ÷ length(x0_range)] + sin(t))
end

# numberOfParameters = length(u0_range)
# function parameterChange(prob, i, repeat)
#    remake(prob, p = t -> u0_range[i] + sin(t))
# end
  
ensembleProb = EnsembleProblem(prob, prob_func=parameterChange)
  
data = solve(ensembleProb, Tsit5(), trajectories=numberOfParameters)
  
using DifferentialEquations.EnsembleAnalysis
#EnsembleSummary(data)
#plot(data)

#temp = transpose(reduce(hcat, data[501].u)) # converts vector of vector to matrix
#scatter(temp[:, 1], temp[:, 2], temp[:, 3]) # can't pass as matrix to scatter and get desired result seemingly

dt = 0.1
t_range = range(t_start, t_final, step=dt)
A = zeros(numberOfParameters, 3, length(t_range))
for (i, t) in enumerate(t_range)
    A[:, :, i] = reduce(hcat, componentwise_vectors_timepoint(data, t))
end

using Interact, Plots, Mux

anim = @gif for (i, t) in enumerate(t_range)
    plot(A[:, 1, i], A[:, 2, i], A[:, 3, i], xlabel = "x", ylabel = "y", zlabel = "cost", title = string(t))
end

plotlyjs()

mp = @manipulate for t=t_range
    i = Int((t - t_start) / dt + 1)
    #scatter(temp[1], temp[2], temp[3], xlabel = "x", ylabel = "y", zlabel = "cost", markersize=0.3)
    plot(A[:, 1, i], A[:, 2, i], A[:, 3, i], xlabel = "x", ylabel = "y", zlabel = "cost")
end
ui = dom"div"(mp);
WebIO.webio_serve(page("/", req -> ui), 8000);

##
println("hello")
##