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
    w2s = collect(range(3702, 4298, length=100)) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.75-4)GHz")
    wcs = collect(range(3000, 4000, length=1000 * 10)) ./ param[:J12].val
    param[:wc] = ValueInfo(wcs, "(3-4.5)GHz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

# function f()
#     param = initial_parameters()
#     # -----------
#     # -----------

#     lines = []
#     d = param[:d]; L = param[:L]

#     zzs = []
#     comparison_states = []

#     println("Calculating zz coupling...")
#     zzs = [zz_zero_point(d, param[:w1].val, w2, param[:wc].val,
#         [param[:U1].val, param[:Uc].val, param[:U2].val],
#         [param[:Jc].val, param[:Jc].val, param[:J12].val]) for w2 in param[:w2].val]
#     zz_min = [minimum(zz) for zz in zzs]
#     zz_max = [maximum(zz) for zz in zzs]
#     g0s = [w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false) for w2 in param[:w2].val]
#     push!(lines, LineInfo(param[:w2].val .* 5, zz_min .* 5, 1, "zz min"))
#     push!(lines, LineInfo(param[:w2].val .* 5, zz_max .* 5, 1, "zz max"))
#     push!(lines, LineInfo(param[:w2].val .* 5, g0s .* 5, 1, "g0"))
    
#     path = joinpath(@__DIR__)
#     makeplot(path, lines...; xlabel = "w2 frequency", ylabel = "wc frequency", copy_code = false, param...)
# end

function add_g0_2()
    param = initial_parameters()
    w2s = param[:w2].val
    g0s = [w_off_position(param[:w1].val - param[:U1].val, w2 - param[:U2].val, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false) for w2 in w2s]
    data = loaddata(joinpath(@__DIR__))

    g0_line = data.lines["g0"]
    zz_min_line = data.lines["zz min"]
    zz_max_line = data.lines["zz max"]
    g0_2_line = LineInfo(w2s .* 5, g0s .* 5, 1, "g0 for leakage")
    pl = plot()
    plot!(pl, g0_line.x, g0_line.y, label = g0_line.tag)
    plot!(pl, zz_min_line.x, zz_min_line.y, label = zz_min_line.tag)
    plot!(pl, zz_max_line.x, zz_max_line.y, label = zz_max_line.tag)
    plot!(pl, g0_2_line.x, g0_2_line.y, label = g0_2_line.tag)
end

# function test_zz_0_interpolation()
#     data = loaddata(joinpath(@__DIR__))

#     zz_min_line = data.lines["zz min"]
#     pl = plot(size=(1000, 800))  # Increase the size of the plot
#     plot!(pl, zz_min_line.x, zz_min_line.y, label = zz_min_line.tag)
#     plot!(pl, 3702:1:4298, x -> zz_0_with_w_4GHZ(x), label = "interpolated", linestyle = :dash)
# end

add_g0_2()

# @time f()