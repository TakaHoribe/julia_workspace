
# import Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("Plots")
# Pkg.add("PlotThemes")
# Pkg.add("PyPlot")

using DifferentialEquations
using LinearAlgebra
using Plots

tau = 0.3

A = [0 1; 0 -1/tau]
B = [0; 1/tau]
u(t) = t < 2.0 ? 1.0 : 0.0

# function exactSolution(x0, tspan, dt)
#     T = tspan[1]:dt:tspan[2]
#     X = [u + exp(-t/tau) * (x0 - u) for t=T]
#     return (T, X)
# end

function calcTsit5(x0, tspan)
    function f(dx, x, p, t)
        dx[1] = x[2]
        dx[2] = -1/tau * x[2] + 1/tau*u(t)
    end

    prob = ODEProblem(f, x0, tspan, [0])
    sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10)
    return (sol.t, mapreduce(permutedims, vcat, sol.u))
end

function integrate(f, x0, tspan, dt)
    x = x0
    T = tspan[1]:dt:tspan[2]
    X = zeros(size(T)[1], 2)

    for i=1:size(T)[1]
        X[i,:] = x
        x = f(x, T[i])
    end
    return (T, X)
end


function taylor1(x0, tspan, dt)
    function f(x, t)
        Ad = I + A*dt
        Bd = B*dt
        x_next = Ad * x + Bd * u(t)
        return x_next
    end

    integrate(f, x0, tspan, dt)
end

function taylor2(x0, tspan, dt)
    function f(x, t)
        Ad = I + A*dt + 1/2*(A*dt)^2
        Bd = (I + 1/2*A*dt) * B * dt
        x_next = Ad * x + Bd * u(t)
        return x_next
    end

    integrate(f, x0, tspan, dt)
end

function bilinear(x0, tspan, dt)
    function f(x, t)
        Ad = inv(I - A*dt/2) * (I + A*dt/2)
        Bd = inv(I - A*dt/2) * B*dt
        x_next = Ad * x + Bd * u(t)
        return x_next
    end

    integrate(f, x0, tspan, dt)
end


x0 = [0.0; 0.0]
tspan = (0.0, 3.0)
dt = 0.2
# sol_exact = exactSolution(x0, tspan, 0.01)
sol_EX = calcTsit5(x0, tspan)
sol_T1 = taylor1(x0, tspan, dt)
sol_T2 = taylor2(x0, tspan, dt)
sol_BL = bilinear(x0, tspan, dt)

lwt = 2
lsw = 3
idx = 2
plot(sol_EX[1], sol_EX[2][:,idx], xaxis="t [s]", yaxis="x", title="\ndiscretization study (closed)", label="Exact Solution", lsw=lsw, linewidth=lwt)
plot!(sol_T1[1], sol_T1[2][:,idx], label="Taylor1", lsw=lsw, linewidth=lwt, linestyle=:dash)
plot!(sol_T2[1], sol_T2[2][:,idx], label="Taylor2", lsw=lsw, linewidth=lwt, linestyle=:dashdot)
plot!(sol_BL[1], sol_BL[2][:,idx], label="Bilinear", lsw=lsw, linewidth=lwt, linestyle=:dashdotdot)
