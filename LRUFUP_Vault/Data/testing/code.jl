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
    param[:d] = 3
    param[:L] = 2
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "BH"
    param[:J] = ValueInfo(0.1, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:U] = ValueInfo(300 / param[:J].val, "300MHz")
    param[:w2] = ValueInfo(4000 / param[:J].val, "4GHz")
    t = 5
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
    t_list = LinRange(0, param[:t].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    e_ops = [n1]

    w = [param[:w1].val, param[:w2].val]
    U = [param[:U].val, param[:U].val]
    H = qt_bosehubbard(d, L; w = w, U = U, J = param[:J].val)
    Ht = H

    lines = []
    state = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0)))
    out = mesolve(Ht, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    line = LineInfo(t_list, real(out.expect[1, :]), 1, "|10>, n1")
    plot(line.x, line.y)
end

@time f()