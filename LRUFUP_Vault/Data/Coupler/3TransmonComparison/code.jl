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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.001
    param[:t] = 50
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(productstate(param[:d], [2, 0, 0]), "|100>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH"
    param[:J12] = ValueInfo(4, "4MHz")
    param[:Jc] = ValueInfo(4 / param[:J12].val, "4MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    #param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:U] = ValueInfo(250 / param[:J12].val, "250MHz")
    param[:U2] = ValueInfo(150 / param[:J12].val, "150MHz")
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

    
    states = [productstate(d, [1, 0, 0]), productstate(d, [2, 0, 0])]
    for state in states
        dU = [param[:U].val, (param[:U].val + param[:U2].val) / 2, param[:U2].val]
        w2 = (2 * param[:w1].val - dU[1] + dU[3]) / 2
        wc = (2 * param[:w1].val - dU[1] + dU[2]) / 2
        ws = [param[:w1].val,  wc, w2]
        H = bosehubbard(d, L, ; w = ws, U = dU)
        Uop = exp(-1im * param[:dt] * Matrix(H))
            
        @time out = exactevolve(state, Uop, param[:dt], param[:t], obs)
        y = out[1, :]
        if state == states[1]
            legend = "Initial state |100>"
        elseif state == states[2]
            legend = "Initial state |200>"
        end
        push!(lines, LineInfo(0:dt:t, y, 1, legend))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "nc", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f_initial_state()