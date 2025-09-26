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
    param[:L] = 2
    N = 2
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)0>")
    param[:alg] = "toolbox"
    param[:H] = "snail"
    param[:phi_c] = ValueInfo(-π / 2, "-π/2")
    param[:Ej] = ValueInfo(20, "20GHz")
    param[:Ec] = ValueInfo(0.3 / param[:Ej].val, "0.3GHz")
    param[:Ej_snail] = ValueInfo(000 / param[:Ej].val, "200GHz")
    param[:N_junc] = 3
    t = 100
    param[:t] = ValueInfo(t * 2 * π * 10^-3 * param[:Ej].val, "$(t)ns")
    param[:diss_rate] = 0

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
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    E_c = [param[:Ec].val, param[:Ec].val]
    E_j = [param[:Ej].val, param[:Ej].val, param[:Ej_snail].val]
    H = qt_two_transmon_snail_test(d, param[:phi_c].val, E_c, E_j, param[:N_junc])

    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n_all"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()