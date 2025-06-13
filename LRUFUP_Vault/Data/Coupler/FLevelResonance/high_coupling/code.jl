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
    param[:state] = ValueInfo(productstate(param[:d], [2, 0, 0]), "|100>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
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

    op = singlesite_n(d, L, 3)
    obs = (state) -> expval(state, op)

    states = [productstate(d, [1, 0, 0]), productstate(d, [2, 0, 0])]

    coupler_coefficient = 100 - 1
    dw = [param[:w].val, param[:w].val + coupler_coefficient, param[:w].val]
    dU = [param[:U].val, param[:U].val + 2 * coupler_coefficient, param[:U].val]
    H = bosehubbard(d, L, [(1, 2, param[:Jc].val), (2, 3, param[:Jc].val), (1, 3, 1)]; w = dw, U = dU)
    Uop = exp(-1im * param[:dt] * Matrix(H))

    for state in states
        @time out = exactevolve(state, Uop, param[:dt], param[:t], obs)
        y = out[1, :]
        if state == states[1]
            legend = "Initial state |100>"
        elseif state == states[2]
            legend = "Initial state |200>"
        end
        push!(lines, LineInfo(0:dt:t, y, 1, legend))
    end

    param[:Uc] = "$(dU[2] * param[:J12].val)MHz"
    excitation_frequencies_plot(dw, dU)

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n2", copy_code = false, param...)
    # addtrajectories(path, line)

end

@time f_initial_state()