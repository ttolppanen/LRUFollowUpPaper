using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays
using Gates

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.01
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(productstate(param[:d], [0, 0, 0]), "|000>")
    param[:alg] = "krylov"
    param[:k] = 6
    param[:H] = "BH"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:U] = ValueInfo(250 / param[:J].val, "Mean 250MHz")
    t = 20
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")
    param[:gate_time] = ValueInfo(10 * 10^(-3) * 2 * π * param[:J].val, "10ns")
    param[:gate_area] = ValueInfo(π, "π")
    param[:gate_repeat_time] = ValueInfo(1 * 2 * π * param[:J].val, "1μs")
    param[:diss_rate] = 0.04
    param[:traj] = 1000
    return param
end

function f()
    param = initial_parameters()

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t].val; k = param[:k]
    wd = param[:w1].val # drive frequency
    state = param[:state].val

    anh_op = anhall(d, L)
    anh_1_op = anhall(d, L; sites = [1])
    obs = [(state) -> expval(state, anh_op),
           (state) -> expval(state, anh_1_op)]


    diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]

    x_op = singlesite_a(d, L, 1) + singlesite_adag(d, L, 1)

    gt = param[:gate_time].val
    ga = param[:gate_area].val

    function time_evolution()
        U = (0.4 * (2 * rand(3) .-1) .+ 1) * param[:U].val
        U = [param[:U].val, U[2], U[3]]
        w2 = (2 * param[:w1].val - U[1] + U[2]) / 2
        w3 = (2 * param[:w1].val - U[1] + U[3]) / 2
        w = [param[:w1].val - wd, w2 - wd, w3 - wd]

        # H = bosehubbard(d, L; w, U, J = 1)
        H = bosehubbard(d, L; w, U, J = 1) .+ H_dd

        effect! = drive_hamiltonian(wd, gt, ga, x_op, H, dt; repeat_time = param[:gate_repeat_time].val, use_rwa = true)
        krylovevolve(state, H, diss_op, dt, t, k, obs...; effect!) 
    end
    
    @time out = solvetrajectories(time_evolution, param[:traj]; paral = :threads)

    out = traj_mean(out)
    push!(lines, LineInfo(0:dt:t, out[1, :], 1, "anh"))
    push!(lines, LineInfo(0:dt:t, out[2, :], 1, "anh_1"))
    # push!(lines, LineInfo(0:dt:t, out[3, :], 1, "anh_3"))
    # push!(lines, LineInfo(0:dt:t, out[4, :], 1, "n1"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()