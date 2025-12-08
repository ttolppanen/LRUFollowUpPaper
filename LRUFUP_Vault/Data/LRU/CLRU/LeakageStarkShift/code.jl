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
    N1 = 0
    Nc = 0
    N2 = 2
    # N3 = 0
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], Nc), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)$(Nc)$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:w2] = ValueInfo(3700 / param[:J].val, "3.7GHz")
    # param[:w3] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:wc1] = ValueInfo(3206 / param[:J].val, "3.206GHz")
    # param[:wc2] = ValueInfo(3206 / param[:J].val, "3.206GHz")
    param[:U] = ValueInfo(300 / param[:J].val, "300MHz")
    t = 0.03
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")

    param[:wd] = ValueInfo(4007 / param[:J].val, "4.007GHz")
    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J].val, "14ns")

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
    anh_1 = qt_single_site_op(d, L, qt_anh_op(d), 1)
    # n3 = qt_single_site_op(d, L, num(d), 5)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2, anh_1]

    # w = [param[:w1].val, param[:wc1].val, param[:w2].val, param[:wc2].val, param[:w3].val]
    w = [param[:w1].val, param[:wc1].val, param[:w2].val] .- param[:wd].val
    U = [param[:U].val for _ in 1:L]
    # g = [param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:J].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J].val]

    # H = qt_tunnel_tunable_coupler(w, U, g)

    I_op, Q_op = qt_drive_op(d, L; rwa = true)
    pulse = gaussian_pulse(param[:gate_time].val, π)

    H = qt_tunable_coupler(w, U, g)
    # H += qt_x_gate_to_QobjEvo(0.0, I_op, Q_op, pulse)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)

    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "anh_1"))
    # push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n2"))
    # push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "n3"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()