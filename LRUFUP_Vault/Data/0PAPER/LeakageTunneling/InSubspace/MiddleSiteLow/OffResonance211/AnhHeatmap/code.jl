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
    ws = collect(range(3700, 3998, 10)) ./ param[:J12].val
    param[:w] = ValueInfo(ws, "(3.7-3.998)GHz")
    param[:w2] = ValueInfo(3700 / param[:J12].val, "3.7GHz")
    param[:wc_min] = ValueInfo(2000 / param[:J12].val, "2GHz")
    param[:wc_max] = ValueInfo(3500 / param[:J12].val, "3.5GHz")
    Us = collect(range(200, 300, 10)) ./ param[:J12].val
    param[:U] = ValueInfo(Us, "200-300MHz")
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

    state211 = productstate(d, [2, 0, 1, 0, 1]); id = find_subspace(state211, indices); state211 = project_to_subspace(state211)
    
    n3 = singlesite_n(d, L, 5); n3 = project_to_subspace(n3)
    e_ops = [QuantumObject(n3[id])]

    max_amps = zeros(length(param[:U].val), length(param[:w].val))
    times = zeros(length(param[:U].val), length(param[:w].val))
    Base.Threads.@threads for I in CartesianIndices((length(param[:U].val), length(param[:w].val)))
        i, j = Tuple(I)
        w3 = param[:w].val[j]
        U_val = param[:U].val[i]

        U_val = param[:U].val[i]

        U = [U_val, U_val, U_val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
        wc_below = NaN
        wc_val = NaN

        try
            wcs, wc_vals = zz_0_with_roots(w3, param[:w2].val, param[:wc_min].val, param[:wc_max].val, U, g)
            wc_below = wcs[1]
            wc_val = wc_vals[1]
        catch
            # println("No ZZ zero point found for w = $(w3 * 5)MHz and U = $(U_val * 5)MHz")
            time = 10 * param[:t].val
            max_amplitude = 0
            max_amps[i, j] = max_amplitude
            times[i, j] = time
            continue
        end

        if abs(wc_val) > 1e-10
            # println("ZZ zero point found for w = $(w3 * 5)MHz and U = $(U_val * 5)MHz, but ZZ is not close to zero (ZZ = $(wc_val))")
            time = 10 * param[:t].val
            max_amplitude = 0
            max_amps[i, j] = max_amplitude
            times[i, j] = time
            continue
        end

        w = [w3, wc_below, param[:w2].val, wc_below, w3]
        U = [U_val, U_val, U...]
        g = [param[:Jc].val, param[:Jc].val, g...]
        
        H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)
        H = QuantumObject(H[id])

        state = QuantumObject(state211[id])
        out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), maxiters = 10^7)
        max_amplitude = maximum(real.(out.expect[1, :]))
        max_i = findfirst(x -> x >= max_amplitude, real.(out.expect[1, :]))
        time = out.times[max_i]
        
        max_amps[i, j] = max_amplitude
        times[i, j] = time
    end
    x = param[:w].val .* 5
    y = param[:U].val .* 5
    # t = t_list ./ (2 * π * 5)
    
    @save joinpath(@__DIR__, "max_amps_data.jld2") x y max_amps
    @save joinpath(@__DIR__, "times_data.jld2") x y times
end

function make_heatmap()
    @load joinpath(@__DIR__, "max_amps_data.jld2") x y max_amps 
    heatmap(x, y, max_amps, xlabel = "w1 = w3 (MHz)", ylabel = "U (MHz)", colorbar_title = "<max(n3)>")
    savefig(joinpath(@__DIR__, "heatmap_max_amps.png"))

    @load joinpath(@__DIR__, "times_data.jld2") x y times
    heatmap(x, y, times ./ (2 * π * 5), clim = (0, 10), color = cgrad(:thermal, rev=true),xlabel = "w1 = w3 (MHz)", ylabel = "U (MHz)", colorbar_title = "t (μs)")
    savefig(joinpath(@__DIR__, "heatmap_times.png"))
end


default_threads = BLAS.get_num_threads()
BLAS.set_num_threads(1)

@time f()
make_heatmap()

BLAS.set_num_threads(default_threads)