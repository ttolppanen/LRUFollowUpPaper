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

    param[:J12] = ValueInfo(8, "5MHz")
    param[:Units] = param[:J12].val
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    U = 400
    offset = 10
    w2s = collect(range(4000 + offset, 4000 - offset, length=100)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.75-4)GHz")
    wcs = collect(range(3000, 4000, length=1000)) ./ param[:J12].val
    param[:wc] = ValueInfo(wcs, "(3-4.5)GHz")
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
    zzs = [zz_zero_point(param[:w1].val, w2, param[:wc].val,
        [param[:U1].val, param[:Uc].val, param[:U2].val],
        [param[:Jc].val, param[:Jc].val, param[:J12].val]) for w2 in param[:w2].val]
    zz_min = [minimum(zz) for zz in zzs]
    zz_max = [maximum(zz) for zz in zzs]
    u = param[:Units]
    g0s = [w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false) for w2 in param[:w2].val]
    push!(lines, LineInfo(param[:w2].val .* u, zz_min .* u, 1, "zz min"))
    push!(lines, LineInfo(param[:w2].val .* u, zz_max .* u, 1, "zz max"))
    push!(lines, LineInfo(param[:w2].val .* u, g0s .* u, 1, "g0"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "w2 frequency", ylabel = "wc frequency", copy_code = false, param...)
end

function plot_difference()
    data = loaddata(joinpath(@__DIR__))

    g0_line = data.lines["g0"]
    zz_min_line = data.lines["zz min"]
    zz_max_line = data.lines["zz max"]
    pl = plot(ylims = [300, 1000])
    plot!(pl, g0_line.x, g0_line.x .- g0_line.y, label = g0_line.tag)
    plot!(pl, zz_min_line.x, zz_min_line.x .- zz_min_line.y, label = zz_min_line.tag)
    plot!(pl, zz_max_line.x, zz_max_line.x .- zz_max_line.y, label = zz_max_line.tag)
end

@time f()
plot_difference()