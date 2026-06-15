using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox
using QuantumTimes
using OrdinaryDiffEq

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 10000
    param[:d] = 4
    param[:L] = 3
    N1 = 1
    N2 = 0
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)0$(N2)>")
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc_below] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:wc_above] = ValueInfo(4495 / param[:J12].val, "4.495GHz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t = 0.1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    # param[:diss_rate] = 0.0

    # param[:wd] = ValueInfo(4007 / param[:J12].val, "4.007GHz")
    # param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J12].val, "14ns")
    # param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:J12].val, "4ns")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    state = param[:state].val

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    # n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1 + nc, n1, qt_n_all(d, L)]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc_below].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]


    H_below = qt_tunable_coupler_paper(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state201 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 1)))
    E, psi, Uni = eigenstates(H_below)
    state10 = psi[4]

    state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))

    state = state10
    out = mesolve(H_below, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |100>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |10>, <N>"))

    state = state201
    out = mesolve(H_below, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |201>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |201>, <N>"))
    
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

function plot_with_analytical()
    plotinfo = loaddata(@__DIR__)
    numerical_line = plotinfo.lines["below, s_0 = |201>, <n1>"]
    pl = plot(numerical_line.x, numerical_line.y, label = numerical_line.tag)

    param = initial_parameters()

    g1c   = param[:Jc].val
    g2c   = -param[:Jc].val          # change sign here if needed by convention
    g12   = param[:J12].val
    Δ     = param[:w1].val - param[:wc_below].val
    U     = param[:U1].val

    gamma2 = 2*g12 + 2*(g1c*g2c)/(Δ - U)

    gamma4 = g12*(g1c^2 + g2c^2)*(-1/Δ^2 - 2/(Δ - U)^2 + 1/(Δ*(Δ - U))) + (g1c*g2c)*(g1c^2 + g2c^2)*(-4/(Δ - U)^3 - 1/(Δ*(Δ - U)^2) + 1/(2*Δ^2*(Δ - U)) + 3/((Δ + U)*(Δ - U)^2) + 1/(Δ*(Δ - U)*(Δ + U)) - 1/((Δ + U)^2*(Δ - U)) - 1/(2*Δ^2*(Δ + U)))

    display(gamma2)
    display(gamma4)
    gamma = gamma2 + gamma4

    plot!(numerical_line.x, x -> cos(gamma * x)^2 + 1, label = "analytical")
    plot!(numerical_line.x, x -> cos(gamma2 * x)^2 + 1, label = "analytical only gamma2")
    plot!(title = PlotAndSave.parameters_to_title(; param..., gamma2, gamma4), titlefont = 8)
    savefig(pl, joinpath(@__DIR__, "analytical_comparison.png"))
end

@time f()
# plot_with_analytical()

