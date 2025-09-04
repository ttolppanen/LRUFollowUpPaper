using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays
using LocalPackages

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function coupler_H(param)
    if haskey(param, :wd)
        return coupler_H(param[:w].val, param[:wd].val, param[:U].val, param[:Jc].val, param[:J12].val, param[:wc].val, param[:Uc].val)
    end
    return coupler_H(param[:w].val, param[:U].val, param[:Jc].val, param[:J12].val, param[:wc].val, param[:Uc].val)
end
function coupler_H(w, U, g_c, g_direct, w_c, U_c)
    coupler_coefficient = g_c^2 / g_direct - 1 # I dont know why there is -1
    # dw = [w, w + coupler_coefficient, w]
    # dU = [U, U + 2.05 * coupler_coefficient, U]
    dw = [w, w_c, w]
    dU = [U, -U_c, U]
    return bosehubbard(3, 3, [(1, 2, -g_c), (2, 3, g_c), (1, 3, g_direct)]; w = dw, U = dU)
end
function coupler_H(w, w_drive, U, g_c, g_direct, w_c, U_c)
    coupler_coefficient = g_c^2 / g_direct - 1 # I dont know why there is -1
    # dw = [w, w + coupler_coefficient, w]
    # dU = [U, U + 2.05 * coupler_coefficient, U]
    dw = [w, w_c, w]
    dU = [U, -U_c, U]
    dw .-= w_drive
    return bosehubbard(3, 3, [(1, 2, -g_c), (2, 3, g_c), (1, 3, g_direct)]; w = dw, U = dU)
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.01
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(productstate(param[:d], [0, 0, 0]), "|000>")
    param[:alg] = "krylov"
    param[:k] = 6
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wd] = ValueInfo(4007 / param[:J12].val, "4.02GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U] = ValueInfo(250 / param[:J12].val, "250MHz")
    param[:Uc] = ValueInfo(732 / param[:J12].val, "740MHz")
    t = 10
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:gate_time] = ValueInfo(20 * 10^(-3) * 2 * π * param[:J12].val, "20ns")
    param[:gate_area] = ValueInfo(π, "π")
    param[:gate_repeat_time] = ValueInfo(0.00 * 2 * π * param[:J12].val, "1μs")
    param[:diss_rate] = 0.4

    param[:J12] = ValueInfo(1, "5MHz")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t].val; k = param[:k]
    state = param[:state].val

    anh_op = anhall(d, L)
    anh_1_op = anhall(d, L; sites = [1])
    anh_3_op = anhall(d, L; sites = [3])
    n1_op = singlesite_n(d, L, 1)
    obs = [(state) -> expval(state, anh_op),
           (state) -> expval(state, anh_1_op),
           (state) -> expval(state, anh_3_op),
           (state) -> expval(state, n1_op)]

    diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]

    # H = coupler_H(param)
    H = coupler_H(param) .+ H_dd

    x_op = singlesite_a(d, L, 1) + singlesite_adag(d, L, 1)
    y_op = 1im * (-singlesite_a(d, L, 1) + singlesite_adag(d, L, 1))

    wd = param[:wd].val
    gt = param[:gate_time].val
    ga = param[:gate_area].val
    effect! = rand_drive_hamiltonian(wd, gt, ga, x_op, y_op, H, dt; repeat_time = param[:gate_repeat_time].val, use_rwa = true)
    
    @time out = krylovevolve(state, H, diss_op, dt, t, k, obs...; effect!)
    push!(lines, LineInfo(0:dt:t, out[1, :], 1, "anh"))
    push!(lines, LineInfo(0:dt:t, out[2, :], 1, "anh_1"))
    push!(lines, LineInfo(0:dt:t, out[3, :], 1, "anh_3"))
    push!(lines, LineInfo(0:dt:t, out[4, :], 1, "n1"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

function vec_to_string(v)
    n1 = singlesite_n(3, 3, 1)
    n2 = singlesite_n(3, 3, 2)
    n3 = singlesite_n(3, 3, 3)
    n1, n2, n3 = expval(v, n1), expval(v, n2), expval(v, n3)
    return "|$(round(n1, digits = 2)) $(round(n2, digits = 2)) $(round(n3))>"
end

function H_ij(H, vecs, i, j)
    println("H_$(i) $(j)")
    v_i = vecs[:, i]
    v_j = vecs[:, j]
    println(abs2(v_i' * H * v_j))
end

function eigen_vals()
    param = initial_parameters()
    H = coupler_H(param)
    E = eigen(Matrix(H))
    pl = scatter(E.values)
    # for i in 1:length(E.values)
    #     println(vec_to_string(E.vectors[:, i]), ": ", E.values[i], "i = $i")
    # end
    H_ij(H, E.vectors, 25, 24)
    pl
end

# @time f()
eigen_vals()