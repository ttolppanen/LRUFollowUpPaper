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
    param[:t] = 0.4
    param[:d] = 2
    param[:L] = 3
    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = 2), "|100...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:J12] = ValueInfo(4, "4MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4100 / param[:J12].val, "4.1GHz")
    param[:w2] = ValueInfo(3890 / param[:J12].val, "3.89GHz")
    param[:U] = ValueInfo(250 / param[:J12].val, "250MHz")
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

    op = singlesite_n(d, L, 1)
    obs = (state) -> expval(state, op)

    for w_coupler in range(2500 / param[:J12].val, 4000 / param[:J12].val, length = 10) 
        ws = [param[:w1].val,  w_coupler, param[:w2].val]
        H = bosehubbard(d, L, [(1, 2, param[:Jc].val), (2, 3, param[:Jc].val), (1, 3, param[:J12].val)]; w = ws, U = [param[:U].val, param[:U].val, -param[:U].val])
        Uop = exp(-1im * param[:dt] * Matrix(H))
            
        @time out = exactevolve(state, Uop, param[:dt], param[:t], obs)
        y = out[1, :]
        push!(lines, LineInfo(0:dt:t, y, 1, "w_coupler = $(w_coupler * param[:J12].val) MHz"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()