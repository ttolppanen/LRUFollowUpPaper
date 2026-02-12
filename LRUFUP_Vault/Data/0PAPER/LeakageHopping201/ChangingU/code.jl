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
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:alg] = "toolbox Vern9"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")

    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    U2s = collect(range(250, 350, length=150)) ./ param[:J12].val
    param[:U2] = ValueInfo(U2s, "(250-350)MHz")

    t = 0.1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
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
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3) - tensor(qeye(d), qeye(d), qeye(d))
    e_ops = [n1 + nc, n1, n2]

    state201 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 1)))

    n2s = []
    for U2 in param[:U2].val
        w2 = (2 * param[:w1].val - param[:U1].val + U2) / 2
        wc = w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false)
        w = [param[:w1].val, wc, w2]
        U = [param[:U1].val, param[:Uc].val, U2]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

        H = qt_tunable_coupler_paper(d, w, U, g)

        state = state201
        out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9())
        push!(n2s, real(out.expect[3, :]))
    end
    n2s = hcat(n2s...)
    x = param[:U2].val .* 5
    t = t_list ./ (2 * π * 5)
    
    heatmap(x, t, n2s, xlabel = "U2 (MHz)", ylabel = "t (μs)", colorbar_title = "<n2-I>")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end


@time f()