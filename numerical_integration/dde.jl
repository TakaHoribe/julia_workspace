using DifferentialEquations
using Plots
gr()

h(p, t) = zeros(2)


function f(dx, x, h, p, t)

    k = p[1]
    b = p[2]
    d = p[3]

    xd = h(p, t-d)
    input = -k * xd[1] - b * xd[2]; # feedback with delay
    p[4] = xd[1]

    dx[1] = x[2]
    dx[2] = input
end

function main(k, b)
    start = time()

    tau = 0.2

    x0 = [1.0, 0.0]
    tspan = (0.0, 30.0)
    param = [k, b, tau, 0.0]

    prob = DDEProblem(f, x0, h, tspan, param);

    saved_values = SavedValues(Float64, Tuple{Float64,Float64});
    cb = SavingCallback((u, t, integrator)->(u[1], integrator.p[4]), saved_values);
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, callback=cb);

    print(saved_values.saveval)
    println("elapsed4 = ", time() - start)

    # popdisplay() # generate new window for plot

    # plot(sol, xaxis="t [s]", yaxis="x", title="ode test", label=["x_1" "x_2"], lsw=3)

    return (sol, saved_values)

end

start = time()
k = 3
b = 1.5
s = SavedValues(Float64, Tuple{Float64,Float64});
(sol, s) = main(k, b);
control_input = [v[2] for v in s.saveval]
plot(sol, xaxis="t [s]", yaxis="x", title="ode test", label=["x_1" "x_2"], lsw=3)
println("elapsed all = ", time() - start);

plot!(sol.t, control_input)