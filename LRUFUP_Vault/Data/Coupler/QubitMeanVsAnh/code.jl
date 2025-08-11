using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function excitation_frequencies_plot(dw, dU)
    first = dw
    second = 2 * dw .- dU
    sites = [1, 2, 3]
    pl = bar(sites, second, label = "Second excitation", xlabel = "Site", ylabel = "Frequency (J12)")
    bar!(pl, sites, first, label = "First excitation")
    savefig(pl, joinpath(@__DIR__, "excitation_frequencies.png"))
    return pl
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.001
    param[:t] = 50
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(productstate(param[:d], [1, 0, 0]), "|100>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:U] = ValueInfo(150 / param[:J12].val, "150MHz")
    return param
end


function f_initial_state()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]
    state = param[:state].val

    op = singlesite_n(d, L, 1)
    obs = (state) -> expval(state, op)

    gcs = collect(5:1:40) ./ param[:J12].val
    
    boson_average = []
    delta_U = []
    for gc in gcs
        coupler_coefficient = gc^2 - 1
        dw = [param[:w].val, param[:w].val + coupler_coefficient, param[:w].val]
        dU = [param[:U].val, param[:U].val + 2 * coupler_coefficient, param[:U].val]
        H = bosehubbard(d, L, [(1, 2, gc), (2, 3, gc), (1, 3, 1)]; w = dw, U = dU)
        Uop = exp(-1im * param[:dt] * Matrix(H))

        out = exactevolve(state, Uop, param[:dt], param[:t], obs)
        push!(boson_average, sum(out[1, :]) / (t / dt))
        push!(delta_U, 2 * coupler_coefficient * param[:J12].val)
    end

    lines = [LineInfo(delta_U, boson_average, 1, "Boson average")]
    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "Delta U MHz", ylabel = "mean population", copy_code = false, param...)
    # addtrajectories(path, line)

    pl = plot(delta_U, boson_average, label = "Boson average", xlabel = "Delta U MHz", ylabel = "mean population")
    plot!(pl, ylims = (0.5, 1))
    hline!(pl, [0.95], label = "0.95")
    savefig(pl, joinpath(@__DIR__, "mean_population_vs_delta_U.png"))
end

@time f_initial_state()