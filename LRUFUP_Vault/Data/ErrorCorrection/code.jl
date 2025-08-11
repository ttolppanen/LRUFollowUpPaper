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
    param[:dt] = 0.02
    param[:t] = 10
    param[:d] = 3
    param[:L] = 3
    param[:alg] = "exact"
    param[:H] = "BH"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")
    return param
end

function f_211()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]
    state = productstate(d, [2, 1, 1])

    s211 = productstate(d, [2, 1, 1])'
    s121 = productstate(d, [1, 2, 1])'
    s112 = productstate(d, [1, 1, 2])'
    obs = [state -> real(s211 * state),
           state -> real(s121 * state),
           state -> real(s112 * state)]

    H = bosehubbard(d, L, ; w = param[:w].val, U = param[:U].val)
    Uop = exp(-1im * param[:dt] * Matrix(H))
        
    @time out = exactevolve(state, Uop, param[:dt], param[:t], obs...)
    push!(lines, LineInfo(0:dt:t, out[1, :], 1, "|211>"))
    push!(lines, LineInfo(0:dt:t, out[2, :], 1, "|121>"))
    push!(lines, LineInfo(0:dt:t, out[3, :], 1, "|112>"))
        
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "overlap", copy_code = false, param...)
    # addtrajectories(path, line)
end

function f_200()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]
    state = productstate(param[:d], [2, 0, 0])

    s211 = productstate(d, [2, 0, 0])'
    s121 = productstate(d, [0, 2, 0])'
    s112 = productstate(d, [0, 0, 2])'
    obs = [state -> real(s211 * state),
           state -> real(s121 * state),
           state -> real(s112 * state)]

    H = bosehubbard(d, L, ; w = param[:w].val, U = param[:U].val)
    Uop = exp(-1im * param[:dt] * Matrix(H))
        
    @time out = exactevolve(state, Uop, param[:dt], param[:t], obs...)
    push!(lines, LineInfo(0:dt:t, out[1, :], 1, "|200>"))
    push!(lines, LineInfo(0:dt:t, out[2, :], 1, "|020>"))
    push!(lines, LineInfo(0:dt:t, out[3, :], 1, "|002>"))
        
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "overlap", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f_211()
@time f_200()