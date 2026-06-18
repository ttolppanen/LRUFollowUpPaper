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
    param[:d] = 5
    param[:L] = 5
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(3750 / param[:J12].val, "3.75GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w3] = ValueInfo(3750 / param[:J12].val, "3.75GHz")
    param[:wc1] = ValueInfo(3162.762 / param[:J12].val, "3.162762GHz")
    param[:wc2] = ValueInfo(3162.762 / param[:J12].val, "3.162762GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 8
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    # U = [param[:U].val, param[:U].val, param[:U].val]
    # gs = [param[:Jc].val, param[:Jc].val, 1]
    # wcs = collect(range(2000, 4200, length=1000)) ./ param[:J12].val
    # wcs = zz_zero_point(param[:w1].val, param[:w2].val, wcs, U, gs)
    # display(wcs .* param[:J12].val)

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
    nc1 = qt_single_site_op(d, L, num(d), 2)
    nc2 = qt_single_site_op(d, L, num(d), 4)
    e_ops = [n1, n2, n3, nc1, nc2, qt_n_all(d, L)]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc1].val, param[:w2].val, param[:wc2].val, param[:w3].val]
    U = [param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val]
    g = [param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:J12].val]

    H = qt_tunnel_tunable_coupler_paper(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state200 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0)))
    state201 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 1)))
    state210 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 0)))
    state211 = QuantumObject(tensor(basis(d, 2), basis(d, 0), basis(d, 1), basis(d, 0), basis(d, 1)))
    state111 = QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0)))
    E, psi, Uni = eigenstates(H)
    state10 = psi[4]

    state = state211
    label_s = "211"
    out = mesolve(H, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "s_0 = |$label_s>, <n1>"))
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "s_0 = |$label_s>, <n2>"))
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "s_0 = |$label_s>, <n3>"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "s_0 = |$label_s>, <nc1>"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "s_0 = |$label_s>, <nc2>"))
    # push!(lines, LineInfo(t_list, real(out.expect[6, :]), 1, "s_0 = |$label_s>, <n tot>"))

    # state = state20
    # out = mesolve(H_below, tate, t_list, e_ops = e_ops, maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "below, s_0 = |200>, <n1>"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

# initial_parameters()

@time f()