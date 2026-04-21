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
    param[:L] = 3
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:alg] = "toolbox Vern9"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    w2s = collect(range(3600, 4000, length=150)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.8-4.2)GHz")

    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t_below = 0.05
    param[:t_below] = ValueInfo(t_below * 2 * π * param[:J12].val, "$(t_below)μs")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    t_list_below = LinRange(0, param[:t_below].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3) - tensor(qeye(d), qeye(d), qeye(d))
    e_ops = [n1, nc, n2]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    state201 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))

    n2s = []
    wcs = []
    for w2 in param[:w2].val
        wc = w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false)
        push!(wcs, wc)
        w = [param[:w1].val, wc, w2]
        U = [param[:U1].val, param[:Uc].val, param[:U2].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

        H = qt_tunable_coupler_paper(d, w, U, g)

        state = state201
        out = mesolve(H, state, t_list_below, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9())
        push!(n2s, real(out.expect[2, :]))
    end
    n2s = hcat(n2s...)
    x = wcs .* 5
    t = t_list_below ./ (2 * π * 5)
    
    heatmap(x, t, n2s, xlabel = "wc (MHz)", ylabel = "t (μs)", colorbar_title = "<nc>")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end


@time f()