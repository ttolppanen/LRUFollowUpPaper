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

# this is calculated so that |2> -> |1> is in resonance
function Ec2(Ej1, Ej2, Ec1)
    c = 2 * (sqrt(8 * Ej1 * Ec1) - Ec1) - Ec1
    Ec2 = (sqrt(32 * Ej2 + 4 * c) - 2 * sqrt(8 * Ej2)) / 2
    return Ec2
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
    param[:phi_c] = ValueInfo(-π * 3 / 2, "-3π/2")
    Ej1 = 20
    param[:Ej1] = ValueInfo(Ej1, "$(Ej1)GHz")
    param[:Ej2] = ValueInfo(15 / 2, "15/2*EJ1")
    param[:Ec1] = ValueInfo(0.3 / param[:Ej1].val, "300MHz")
    param[:Ec2] = ValueInfo(Ec2(1, param[:Ej2].val, param[:Ec1].val), "|2>->|1> in resonance")
    param[:Ej_snail] = ValueInfo(3 / 2 , "3/2*EJ1")
    param[:N_junc] = 2
    t = 100
    param[:t] = ValueInfo(t * 2 * π * 10^-3 * param[:Ej1].val, "$(t)ns")
    param[:diss_rate] = 0

    param[:Ej1] = ValueInfo(1, param[:Ej1].str_val)
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

    E_c = [param[:Ec1].val, param[:Ec2].val]
    E_j = [param[:Ej1].val, param[:Ej2].val, param[:Ej_snail].val]
    H = qt_two_transmon_snail(d, param[:phi_c].val, E_c, E_j, param[:N_junc])

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