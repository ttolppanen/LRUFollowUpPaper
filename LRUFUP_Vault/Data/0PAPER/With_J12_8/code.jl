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
    param[:t_points] = 1000
    param[:d] = 4
    param[:L] = 3
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(8, "8MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc_below] = ValueInfo(3694.5 / param[:J12].val, "3.7Hz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:diss_rate] = 0.1

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
    n2 = qt_single_site_op(d, L, num(d), 3)
    n_all = qt_n_all(d, L)
    e_ops = [n1 + nc, n1, nc, n2]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc_below].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]


    H_below = qt_tunable_coupler_paper(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state200 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    state201 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 1)))
    E, psi, Uni = eigenstates(H_below)

    # display(length(psi))
    # for i in 1:length(psi)
    #     s = psi[i]
    #     val = real(expect(n_all, s))
    #     if val < 2
    #         println("hep $i $val")
    #     end
    # end
    # state10 = psi[3]

    # state = state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))
    # out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "below, s_0 = |100>, <n1 + nc>"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |100>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "below, s_0 = |100>, <n2>"))

    # state = state200
    # out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |200>, <n1>"))
    
    # state = state201
    # out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |201>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |201>, <nc>"))

    # state = state010 = QuantumObject(tensor(basis(param[:d], 0), basis(param[:d], 1), basis(param[:d], 0)))
    # out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |010>, <nc>"))

    state = state011 = QuantumObject(tensor(basis(param[:d], 0), basis(param[:d], 1), basis(param[:d], 1)))
    out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "below, s_0 = |011>, <nc>"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |011>, <n1>"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "below, s_0 = |011>, <n2>"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()


