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
    param[:L] = 9
    N1 = 1
    N2 = 1
    N3 = 1
    N4 = 2
    N5 = 1
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2), basis(param[:d], 0), basis(param[:d], N3), basis(param[:d], N4) , basis(param[:d], 0) , basis(param[:d], 0) , basis(param[:d], N5)))
    param[:state] = ValueInfo(state, "|$(N1)$(N2)$(N3)$(N4)$(N5)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J].val, "50MHz")
    param[:w_high] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:w_low] = ValueInfo(3700 / param[:J].val, "3.7GHz")
    param[:wc] = ValueInfo(3200 / param[:J].val, "3.2GHz")
    param[:U] = ValueInfo(300 / param[:J].val, "300MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")

    param[:J] = ValueInfo(1, param[:J].str_val)
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
    n3 = qt_single_site_op(d, L, num(d), 5)
    n4 = qt_single_site_op(d, L, num(d), 6)
    n5 = qt_single_site_op(d, L, num(d), 9)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2, n3, n4, n5]

    wh = param[:w_high].val
    wl = param[:w_low].val
    wc = param[:wc].val
    w = [wl, wc, wh, wc, wl, wl, wc, wc, wl]
    U = [param[:U].val for _ in 1:L]

    H = qt_cross_tunable_coupler(w, U, param[:Jc].val, param[:J].val)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)

    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n2"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "n3"))
    push!(lines, LineInfo(t_list, real(out.expect[6, :]), 1, "n4"))
    push!(lines, LineInfo(t_list, real(out.expect[7, :]), 1, "n5"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()