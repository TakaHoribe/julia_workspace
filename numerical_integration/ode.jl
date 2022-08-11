
# import Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("Plots")
# Pkg.add("PyPlot")

# using DifferentialEquations
# using Plots


k = 1
b = 0.0

function f(dx, x, p, t)
    k = p[1]
    b = p[2]
    dx[1] = x[2]
    dx[2] = -k * x[1] - b * x[2]
end

x0 = [1.0, 1.0]
tspan = (0.0, 30.0)
param = [k, b]

prob = ODEProblem(f, x0, tspan, param)
start = time()
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-3)
println("elapsed1 = ", time() - start)

plot(sol, xaxis="t [s]", yaxis="x", title="ode test", label=["x_1" "x_2"], lsw=3)
