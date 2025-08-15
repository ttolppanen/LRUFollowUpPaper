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

function drive_hamiltonian(H, dt, drive_strength, op, cut_off_time)
    total_time = [0.0]
    H0 = deepcopy(H)

    function Ht(state)
        t = total_time[1]
        if t < cut_off_time
            t += dt
            H .= H0 - drive_strength * cos(t) * op 
        end
        return state
    end
    return Ht
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:t] = 50
    param[:d] = 3
    param[:L] = 1
    s0 = complex(spzeros(param[:d]))
    s0[1] = 1.0
    param[:state] = ValueInfo(s0, "|0>")
    param[:alg] = "krylov"
    param[:k] = 6
    param[:H] = "BH"
    param[:w] = ValueInfo(4000, "4GHz")
    param[:U] = ValueInfo(200 / param[:w].val, "200MHz")

    param[:w] = ValueInfo(1, "4GHz")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]; k = param[:k]
    state = param[:state].val

    op = nop(d)
    obs = (state) -> expval(state, op)

    H = bosehubbard(d, L; w = param[:w].val, U = param[:U].val)

    x_op = aop(d) + adagop(d)
    effect! = drive_hamiltonian(H, dt, 1, x_op, 1)
    @time out = krylovevolve(state, H, dt, t, k, obs; effect!)
    y = out[1, :]
    push!(lines, LineInfo(0:dt:t, y, 1, "with drive"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()