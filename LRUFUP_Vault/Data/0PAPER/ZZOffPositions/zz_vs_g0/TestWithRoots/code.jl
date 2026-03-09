using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox


function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3

    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    w2s = collect(range(3727, 4000, length=100)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.75-4)GHz")
    param[:min_wc] = ValueInfo(3000 / param[:J12].val, "3GHz")
    param[:max_wc] = ValueInfo(4000 / param[:J12].val, "4GHz")
    U = 275
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    zzs = []
    comparison_states = []

    println("Calculating zz coupling...")
    zzs = [zz_0_with_roots(param[:w1].val, w2, param[:min_wc].val, param[:max_wc].val,
        [param[:U1].val, param[:Uc].val, param[:U2].val],
        [param[:Jc].val, param[:Jc].val, param[:J12].val]) for w2 in param[:w2].val]
    zz_min = [zz[1][1] for zz in zzs]
    zz_max = [zz[1][2] for zz in zzs]
    zz_min_values = [zz[2][1] for zz in zzs]
    zz_max_values = [zz[2][2] for zz in zzs]
    # g0s = [w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false) for w2 in param[:w2].val]
    push!(lines, LineInfo(param[:w2].val .* 5, zz_min .* 5, 1, "zz min"))
    # push!(lines, LineInfo(param[:w2].val .* 5, zz_max .* 5, 1, "zz max"))
    # push!(lines, LineInfo(param[:w2].val .* 5, g0s .* 5, 1, "g0"))

    # comparison_states = []
    # zz_min_values = [calc_zz(zz_min[i], comparison_states, param[:w1].val, param[:w2].val[i], [param[:U1].val, param[:Uc].val, param[:U2].val], [param[:Jc].val, param[:Jc].val, param[:J12].val]) for i in eachindex(zz_min)]
    # zz_max_values = [calc_zz(zz_max[i], comparison_states, param[:w1].val, param[:w2].val[i], [param[:U1].val, param[:Uc].val, param[:U2].val], [param[:Jc].val, param[:Jc].val, param[:J12].val]) for i in eachindex(zz_max)]
    pl = plot()
    plot!(pl, param[:w2].val .* 5, zz_min_values; label = "zz at min")
    plot!(pl, param[:w2].val .* 5, zz_max_values; label = "zz at max")
    savefig(pl, joinpath(@__DIR__, "zz_at_extrema.png"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "w2 frequency", ylabel = "wc frequency", copy_code = false, param...)
end

@time f()