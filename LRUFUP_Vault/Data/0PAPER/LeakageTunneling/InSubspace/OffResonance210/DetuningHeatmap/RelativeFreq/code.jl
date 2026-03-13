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
using JLD2

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
    ws = collect(range(3802, 4000, 100)) ./ param[:J12].val
    param[:w] = ValueInfo(ws, "(3.802-4)GHz")
    offset_N = 50
    w_offsets = vcat(collect(range(-50, 0, offset_N)), collect(range(0, 50, offset_N)))
    param[:w_offset] = ValueInfo(w_offsets ./ param[:J12].val, "(-50-50)MHz")
    param[:wc_min] = ValueInfo(2500 / param[:J12].val, "3GHz")
    param[:wc_max] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 10
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

    state210 = productstate(d, [2, 0, 1, 0, 0]); id = find_subspace(state210, indices); state210 = project_to_subspace(state210)
    state = QuantumObject(state210[id])

    n3 = singlesite_n(d, L, 5); n3 = project_to_subspace(n3)
    e_ops = [QuantumObject(n3[id])]

    x_N = length(param[:w].val)
    y_N = length(param[:w_offset].val)
    max_amps = zeros(y_N, x_N)
    Base.Threads.@threads for I in CartesianIndices((y_N, x_N))
        i, j = Tuple(I)
        w1 = param[:w].val[j]
        w3 = param[:w_offset].val[i] + w1

        U = [param[:U].val, param[:U].val, param[:U].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
        wc_below_1c = NaN
        wc_val_1c = NaN
        wc_below_2c = NaN
        wc_val_2c = NaN

        try
            wcs, wc_vals = zz_0_with_roots(w1, param[:w2].val, param[:wc_min].val, param[:wc_max].val, U, g)
            wc_below_1c = wcs[1]
            wc_val_1c = wc_vals[1]

            wcs, wc_vals = zz_0_with_roots(w3, param[:w2].val, param[:wc_min].val, param[:wc_max].val, U, g)
            wc_below_2c = wcs[1]
            wc_val_2c = wc_vals[1]
        catch
            # println("No ZZ zero point found for w = $(w3 * 5)MHz and U = $(U_val * 5)MHz")
            max_amplitude = 0
            max_amps[i, j] = max_amplitude
            continue
        end

        if abs(wc_val_1c) > 1e-10 || abs(wc_val_2c) > 1e-10
            # println("ZZ zero point found for w = $(w3 * 5)MHz and U = $(U_val * 5)MHz, but ZZ is not close to zero (ZZ = $(wc_val))")
            max_amplitude = 0
            max_amps[i, j] = max_amplitude
            continue
        end

        w = [w1, wc_below_1c, param[:w2].val, wc_below_2c, w3]
        U = [param[:U].val, param[:U].val, U...]
        g = [param[:Jc].val, param[:Jc].val, g...]
        
        H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)
        H = QuantumObject(H[id])

        out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), maxiters = 10^7)
        max_amplitude = maximum(real.(out.expect[1, :]))
        max_i = findfirst(x -> x >= max_amplitude, real.(out.expect[1, :]))
        
        max_amps[i, j] = max_amplitude
    end
    x = param[:w].val .* 5
    y = param[:w_offset].val .* 5
    # t = t_list ./ (2 * π * 5)
    
    @save joinpath(@__DIR__, "max_amps_data.jld2") x y max_amps
end

function make_heatmap()
    @load joinpath(@__DIR__, "max_amps_data.jld2") x y max_amps
    heatmap(x, y, max_amps, xlabel = "w1 (MHz)", ylabel = "w3 - w1 (MHz)", colorbar_title = "<max(n3)>")
    savefig(joinpath(@__DIR__, "heatmap_max_amps.png"))
end


@time f()
make_heatmap()