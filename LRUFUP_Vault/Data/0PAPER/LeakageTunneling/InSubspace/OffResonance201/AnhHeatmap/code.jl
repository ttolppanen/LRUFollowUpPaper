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
using Colors

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
    ws = collect(range(3802, 4000, 10)) ./ param[:J12].val
    param[:w] = ValueInfo(ws, "(3.752-4)GHz")
    param[:wc_min] = ValueInfo(3000 / param[:J12].val, "3GHz")
    param[:wc_max] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    Us = collect(range(200, 300, 2)) ./ param[:J12].val
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

    state201 = productstate(d, [2, 0, 0, 0, 1]); id = find_subspace(state201, indices); state201 = project_to_subspace(state201)
    
    n3 = singlesite_n(d, L, 5); n3 = project_to_subspace(n3)
    e_ops = [QuantumObject(n3[id])]

    max_amps_time = []
    for w in param[:w].val
        U_list = []
        for U_val in param[:U].val
            U = [U_val, U_val, U_val]
            g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
            wc_below = NaN
            wc_val = NaN

            try
                wcs, wc_vals = zz_0_with_roots(w, param[:w2].val, param[:wc_min].val, param[:wc_max].val, U, g)
                wc_below = wcs[1]
                wc_val = wc_vals[1]
            catch
                time = 10 * param[:t].val
                max_amplitude = 0
                push!(U_list, (max_amplitude, time))
                continue
            end

            if wc_val > 1e-10
                time = 10 * param[:t].val
                max_amplitude = 0
                push!(U_list, (max_amplitude, time))
                continue
            end

            w = [w, wc_below, param[:w2].val, wc_below, w]
            U = [U_val, U_val, U...]
            g = [param[:Jc].val, param[:Jc].val, g...]
            
            H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)
            H = QuantumObject(H[id])

            state = QuantumObject(state201[id])
            out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), maxiters = 10^7)
            max_amplitude = maximum(real.(out.expect[1, :]))
            max_i = findfirst(x -> x >= max_amplitude, real.(out.expect[1, :]))
            time = out.times[max_i]
            push!(U_list, (max_amplitude, time))
        end
        push!(max_amps_time, U_list)
    end
    max_amps_time = hcat(max_amps_time...)
    x = param[:w].val .* 5
    y = param[:U].val .* 5
    # t = t_list ./ (2 * π * 5)

    cmap = cgrad(:RdYlGn)
    colors = similar(max_amps_time, RGB{Float64})

    for y_i in axes(max_amps_time, 1), x_i in axes(max_amps_time, 2)
        amp = max_amps_time[y_i, x_i][1]
        time = max_amps_time[y_i, x_i][2]
        println(amp)
        println(time)
        c = cmap[time / param[:t].val] 
        norm_amp = amp / 2
        c = RGB(norm_amp * c.r, norm_amp * c.g, norm_amp * c.b)
        colors[y_i, x_i] = c
    end

    display(colors)
    
    heatmap(x, y, color = colors, xlabel = "w1 = w3 (MHz)", ylabel = "U (MHz)", colorbar_title = "<t (μs)>")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end


@time f()