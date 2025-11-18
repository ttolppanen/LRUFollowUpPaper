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

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    N1 = 2
    N2 = 1
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)0$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler +J"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U1] = ValueInfo(300 / param[:J12].val, "300MHz")
    U2 = 300
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    param[:Uc] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 0.1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:diss_rate] = 0.0

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
    n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), 1)]

    w = [param[:w1].val, param[:wc].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    println("J = $(J_eff(g, w))")

    # I_op, Q_op = qt_drive_op(d, L; rwa = true)
    # pulse = gaussian_pulse(param[:gate_time].val, π)

    H = qt_tunable_coupler(w, U, g)
    # for t_offset in [2 * i for i in 0:20]
    #     H += qt_x_gate_to_QobjEvo(0.0, I_op, Q_op, t -> pulse(t - t_offset))
    # end

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))

    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    n2_line = LineInfo(t_list, real(out.expect[4, :]), 1, "n2")
    push!(lines, n2_line)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
    w_second_level = [param[:w1].val - param[:U1].val, param[:wc].val - 1, param[:w2].val - param[:U2].val]
    
    Δ = w_second_level[1] - w_second_level[3]
    pl = plot(n2_line.x, n2_line.y, label = "Δ = $(Δ * 5)")
    J_e = J_eff_new(w_second_level, g) * 1.62
    println("Effective J in second level: ", round(J_e * 5, digits=6), " MHz")
    JJ = 4 * J_e^2
    omega = sqrt(Δ^2 + JJ)
    plot!(pl, n2_line.x, t -> JJ / omega^2 * sin(omega / 2 * t)^2 + 1, label = "Δ = $(Δ * 5), rabi")
    savefig(path * "/plot.png")
end

function t_time_from_n(out)
    x = out.times
    y = real(out.expect[1, :])
    pl = plot(x, y)
    n = length(out.times) # adjust here to make fit better
    c, slope, f_fit = exp_fit(x[1:n], y[1:n], 2)
    println(-1 / slope / (2 * π * 5 * 1e6 * 1e-9), "ns")
    plot!(x, f_fit)
end

function J_eff(g, w)
    delta_1 = w[1] - w[2] + 1
    delta_2 = w[3] - w[2] + 1
    one_over_delta = 0.5 * (1 / delta_1 + 1 / delta_2)
    return -g[1] * g[2] * one_over_delta + g[3]
end

function J_eff_new(w, g)
    delta_1 = w[1] - w[2]
    delta_2 = w[3] - w[2]
    one_over_delta = 0.5 * (1 / delta_1 + 1 / delta_2)
    j_eff = -g[1] * g[2] * one_over_delta + g[3]
    if abs(j_eff) > (g[1] + g[3])
        return (g[1] + g[3])
    end
    return j_eff
end

@time pl = f()