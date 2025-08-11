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
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    #param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:U] = ValueInfo(240 / param[:J12].val, "240MHz")
    param[:U2] = ValueInfo(200 / param[:J12].val, "200MHz")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]
    state = param[:state].val

    op = singlesite_n(d, L, 3)
    obs = (state) -> expval(state, op)

    for w_coupler in [0, 3875, 4625] ./ param[:J12].val
        w2 = (2 * param[:w1].val - param[:U].val + param[:U2].val) / 2
        ws = [param[:w1].val,  w_coupler, w2]
        H = bosehubbard(d, L, [(1, 2, param[:Jc].val), (2, 3, param[:Jc].val), (1, 3, 1)]; w = ws, U = [param[:U].val, param[:U].val, param[:U2].val])
        Uop = exp(-1im * param[:dt] * Matrix(H))
            
        @time out = exactevolve(state, Uop, param[:dt], param[:t], obs)
        y = out[1, :]
        push!(lines, LineInfo(0:dt:t, y, 1, "w_coupler = $(w_coupler * param[:J12].val) MHz"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
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

    w_coupler = 3900 / param[:J12].val
    param[:wc] = "0GHz"
    states = [productstate(d, [1, 0, 0]), productstate(d, [2, 0, 0])]
    for state in states
        w2 = (2 * param[:w1].val - param[:U].val + param[:U2].val) / 2
        ws = [param[:w1].val,  w_coupler, w2]
        H = bosehubbard(d, L, [(1, 2, param[:Jc].val), (2, 3, param[:Jc].val), (1, 3, 1)]; w = ws, U = [param[:U].val, param[:U].val, param[:U2].val])
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
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f_initial_state()