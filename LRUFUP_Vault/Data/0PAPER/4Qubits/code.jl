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
    param[:L] = 7
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(3750 / param[:J12].val, "3.75GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w3] = ValueInfo(4000 / param[:J12].val, "3.75GHz")
    param[:w4] = ValueInfo(4000 / param[:J12].val, "3.75GHz")
    param[:wc1] = ValueInfo(3642 / param[:J12].val, "3.713313GHz")
    param[:wc2] = ValueInfo(3505 / param[:J12].val, "3.713313GHz")
    param[:wc3] = ValueInfo(3505 / param[:J12].val, "3.713313GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 0.1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    U = [param[:U].val, param[:U].val, param[:U].val]
    gs = [param[:Jc].val, param[:Jc].val, param[:J12].val]
    wcs = collect(range(2000, 4200, length=1000)) ./ param[:J12].val
    wcs = zz_zero_point(param[:w1].val, param[:w2].val, wcs, U, gs)
    display(wcs .* param[:J12].val)

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
    n2 = qt_single_site_op(d, L, num(d), 3)
    n3 = qt_single_site_op(d, L, num(d), 5)
    n4 = qt_single_site_op(d, L, num(d), 7)
    nc1 = qt_single_site_op(d, L, num(d), 2)
    nc2 = qt_single_site_op(d, L, num(d), 4)
    nc3 = qt_single_site_op(d, L, num(d), 6)
    e_ops = [n1, n2, n3, n4, nc1, nc2, nc3]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc1].val, param[:w2].val, param[:wc2].val, param[:w3].val, param[:wc3].val, param[:w4].val]
    U = [param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val]
    g = [param[:Jc].val, param[:J12].val]

    H = qt_four_qubit_coupler_paper(w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state2000 = QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 0)))
    state2100 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0)))
    # state210 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0)))
    # state211 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 0)))
    E, psi, Uni = eigenstates(H)
    state10 = psi[4]

    state = state2000
    label_s = "2000"
    out = mesolve(H, state, t_list, e_ops = e_ops, maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "s_0 = |$label_s>, <n1>"))
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "s_0 = |$label_s>, <n2>"))
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "s_0 = |$label_s>, <n3>"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "s_0 = |$label_s>, <n4>"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "s_0 = |$label_s>, <nc1>"))
    push!(lines, LineInfo(t_list, real(out.expect[6, :]), 1, "s_0 = |$label_s>, <nc2>"))
    push!(lines, LineInfo(t_list, real(out.expect[7, :]), 1, "s_0 = |$label_s>, <nc3>"))

    # state = state20
    # out = mesolve(H_below, state, t_list, e_ops = e_ops, maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |200>, <n1>"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

# initial_parameters()

@time f()