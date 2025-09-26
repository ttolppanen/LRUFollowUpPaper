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
    param[:L] = 1
    N = 1
    state = basis(param[:d], N)
    param[:state] = ValueInfo(state, "|$(N)>")
    param[:alg] = "toolbox"
    param[:H] = "snail"
    param[:phi_c] = ValueInfo(-pi / 2, "-π/2")
    param[:Ej] = ValueInfo(20, "20GHz")
    param[:Ec] = ValueInfo(0.3 / param[:Ej].val, "0.3GHz")
    t = 30
    param[:t] = ValueInfo(t * 2 * π * 10^-3 * param[:Ej].val, "$(t)ns")
    # param[:diss_rate] = 1

    param[:Ej] = ValueInfo(1, param[:Ej].str_val)
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
    n_all = qt_n_all(d, L)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n_all]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    E_c = [param[:Ec].val, param[:Ec].val]
    E_j = [param[:Ej].val, param[:Ej].val, param[:Ej].val]
    H = qt_single_transmon_snail(d, param[:phi_c].val, E_c, E_j, 1)
    display(H)

    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n_all"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()