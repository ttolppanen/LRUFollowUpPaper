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
    param[:w1] = ValueInfo(3750 / param[:J12].val, "3.75GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    length = 100
    w3s = vcat(collect(range(3745, 3750, length)), collect(range(3750, 3755, length))) ./ param[:J12].val
    param[:w3] = ValueInfo(w3s, "(3.995-4.005)GHz")
    param[:wc1] = ValueInfo(3162.762 / param[:J12].val, "3.162762GHz")
    param[:wc2] = ValueInfo(3162.762 / param[:J12].val, "3.162762GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 2
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    # U = [param[:U].val, param[:U].val, param[:U].val]
    # gs = [param[:Jc].val, param[:Jc].val, 1]
    # wcs = collect(range(2000, 4200, length=1000)) ./ param[:J12].val
    # wcs = zz_zero_point(param[:w1].val, param[:w2].val, wcs, U, gs)
    # display(wcs .* param[:J12].val)

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

    max_amplitude_below = []
    for w3 in param[:w3].val
        # wc = w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false)
        w = [param[:w1].val, param[:wc1].val, param[:w2].val, param[:wc2].val, w3]
        U = [param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val]
        g = [param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:J12].val]

        H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)
        H = QuantumObject(H[id])

        state = QuantumObject(state211[id])
        out = mesolve(H, state, t_list, e_ops = e_ops, progress_bar= false; abstol = param[:abstol], reltol = param[:reltol], alg = Vern9(), maxiters = 10^7)
        push!(max_amplitude_below, maximum(real(out.expect[1, :])))
    end
    push!(lines, LineInfo(param[:w3].val .* 5, max_amplitude_below, 1, "max(<n3>)"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "MHz", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()