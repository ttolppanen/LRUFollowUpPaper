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
    param[:d] = 5
    param[:L] = 5
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    ws = collect(range(3800, 4000, 200)) ./ param[:J12].val
    param[:w] = ValueInfo(ws, "(3.752-4)GHz")
    wcs = collect(range(3000, 4000, 1000)) ./ param[:J12].val
    param[:wcs] = ValueInfo(wcs, "(3-4)GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:U] = ValueInfo(275 / param[:J12].val, "275MHz")
    t = 8
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

    indices = total_boson_number_subspace_indices(d, L)
    project_to_subspace(m) = subspace_split(m, indices)

    state200 = productstate(d, [2, 0, 0, 0, 0]); id = find_subspace(state200, indices); state200 = project_to_subspace(state200)
    
    n3 = singlesite_n(d, L, 5); n3 = project_to_subspace(n3)
    e_ops = [QuantumObject(n3[id])]

    n3s = []
    for w in param[:w].val
        U = [param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val]
        g = [param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:J12].val]
        wc = minimum(zz_zero_point(w, param[:w2].val, param[:wcs].val, U, g))
        w = [w, wc, param[:w2].val, wc, w]

        H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)
        H = QuantumObject(H[id])

        state = QuantumObject(state200[id])
        out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), maxiters = 10^7)
        push!(n3s, real(out.expect[1, :]))
    end
    n3s = hcat(n3s...)
    x = param[:w].val .* 5
    t = t_list ./ (2 * π * 5)
    
    heatmap(x, t, n3s, xlabel = "w1 = w3 (MHz)", ylabel = "t (μs)", colorbar_title = "<n3>")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end


@time f()