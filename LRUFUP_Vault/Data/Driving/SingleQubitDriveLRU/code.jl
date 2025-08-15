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

function drive_hamiltonian(H, w, dt, drive_strength, op, cut_off_time)
    total_time = [0.0]
    H0 = deepcopy(H)
    println(w)

    function Ht(state)
        t = total_time[1]
        if t < cut_off_time
            t += dt
            H .= H0 - drive_strength / 2 * cos(w * t) * op 
        end
        return state
    end
    return Ht
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.005
    param[:t] = 100
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(productstate(param[:d], [0, 0, 0]), "|000>")
    param[:alg] = "krylov"
    param[:k] = 6
    param[:H] = "BH"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:U] = ValueInfo(250 / param[:J].val, "Mean 250MHz")
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

    op = singlesite_n(d, L, 1)
    obs = (state) -> expval(state, op)

    U = (0.2 * (2 * rand(3) .-1) .+ 1) * param[:U].val
    w2 = (2 * param[:w1].val - U[1] + U[2]) / 2
    w3 = (2 * param[:w1].val - U[1] + U[3]) / 2
    w = [param[:w1].val, w2, w3]

    H = bosehubbard(d, L; w, U, J = 1)

    x_op = singlesite_a(d, L, 1) + singlesite_adag(d, L, 1)
    effect! = drive_hamiltonian(H, param[:w1].val, dt, 1, x_op, 1)
    @time out = krylovevolve(state, H, dt, t, k, obs; effect!)
    y = out[1, :]
    push!(lines, LineInfo(0:dt:t, y, 1, "with drive"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()