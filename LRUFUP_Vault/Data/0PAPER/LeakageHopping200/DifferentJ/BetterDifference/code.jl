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
    param[:d] = 3
    param[:L] = 3
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(45 / param[:J12].val, "45MHz")
    param[:Jc2] = ValueInfo(55 / param[:J12].val, "55MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(3995.6 / param[:J12].val, "3.998365GHz")
    param[:wc_below] = ValueInfo(3495 / param[:J12].val, "3.495GHz")
    # param[:wc_above] = ValueInfo(4495 / param[:J12].val, "4.495GHz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t = 2
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

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    # n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1 + nc, n1, qt_n_all(d, L)]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc_below].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc2].val, param[:J12].val]


    H_below = qt_tunable_coupler_paper(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    E, psi, Uni = eigenstates(H_below)
    state10 = psi[4]

    state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))

    state = state10
    out = mesolve(H_below, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |100>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |10>, <N>"))

    state = state20
    out = mesolve(H_below, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |200>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |200>, <N>"))
    
    # w = [param[:w1].val, param[:wc_above].val, param[:w2].val]
    # H_above = qt_tunable_coupler_paper_positive(w, U, g)

    # E, psi, Uni = eigenstates(H_above)
    # state10 = psi[3]

    # state = state10
    # out = mesolve(H_above, state, t_list, e_ops = e_ops)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "above, s_0 = |10>, <nc+n1>"))

    # state = state20
    # out = mesolve(H_above, state, t_list, e_ops = e_ops, maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "above, s_0 = |200>, <n1>"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

function plot_with_analytical()
    plotinfo = loaddata(@__DIR__)
    numerical_line = plotinfo.lines["below, s_0 = |200>, <n1>"]
    pl = plot(numerical_line.x, numerical_line.y, label = numerical_line.tag)

    param = initial_parameters()
    g1c = param[:Jc].val
    g2c = -param[:Jc].val
    g12 = param[:J12].val
    delta = param[:w1].val - param[:wc_below].val
    U = param[:Uc].val
    gamma1 = 2*g12*g1c*g2c*(1/(delta-U)^2 - 1/(delta*(delta-U))) + 0.5*g1c^2*g2c^2*(-3/(delta*(delta-U)^2) + 6/((delta-U)^2*(delta+U)) - 1/(delta^2*(delta-U)) - 2/((delta-U)*(delta+U)^2))

    gamma2 = 2*(g12 + (g1c * g2c / 2) * (1 / (delta - U) + 1 / delta))^2 / U
    
    display(gamma1)
    display(gamma2)

    plot!(numerical_line.x, x -> cos(x * (gamma1 + gamma2))^2 * 2, label = "analytical")
    plot!(numerical_line.x, x -> cos(gamma1 * x)^2 * 2, label = "analytical only gamma1")
    plot!(title = PlotAndSave.parameters_to_title(; param..., gamma1, gamma2), titlefont = 8)
    savefig(joinpath(@__DIR__, "analytical comparison"))
end

@time f()

# plot_with_analytical()

