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
    N1 = 1
    N2 = 0
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)0$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U1] = ValueInfo(300 / param[:J12].val, "300MHz")
    U2 = 300
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    param[:Uc] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 2
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:diss_rate] = 0.0

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
    n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2, nc, n1 + nc]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    println("J = $(J_eff(g, w))")

    # I_op, Q_op = qt_drive_op(d, L; rwa = true)
    # pulse = gaussian_pulse(param[:gate_time].val, π)

    H = qt_tunable_coupler_below(w, U, g)
    # for t_offset in [2 * i for i in 0:20]
    #     H += qt_x_gate_to_QobjEvo(0.0, I_op, Q_op, t -> pulse(t - t_offset))
    # end

    t_list = LinRange(0, param[:t].val, param[:t_points])
    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    state210 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 1), basis(param[:d], 0)))
    state010 = QuantumObject(tensor(basis(param[:d], 0), basis(param[:d], 1), basis(param[:d], 0)))
    state11 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 1), basis(param[:d], 0)))

    E, psi, U = eigenstates(H)
    state10 = psi[4]

    state = state10
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "initial_state = logical |10>, <n1>"))
    push!(lines, LineInfo(t_list, real(out.expect[6, :]), 1, "initial_state = logical |10>, <nc + n1>"))

    state = state20
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "initial_state = |200>, <n1>"))

    state = state11
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "initial_state = |110>, <n1>"))
     push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "initial_state = |110>, <n2>"))

    # state = state010
    # out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    # push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "initial_state = |010>, <n1>"))

    state = state210
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "initial_state = |210>, <n1>"))

    # out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))

    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n2"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
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

@time f()