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
    param[:wd] = ValueInfo(4009.1 / param[:J12].val, "4.0091GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U] = ValueInfo(250 / param[:J12].val, "250MHz")
    param[:Uc] = ValueInfo(732 / param[:J12].val, "732MHz")
    param[:gate_time] = ValueInfo(20 * 10^(-3) * 2 * π * param[:J12].val, "20ns")

    param[:J12] = ValueInfo(1, "5MHz")
    return param
end

function f()
    param = initial_parameters()
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 1
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    dt = param[:dt]; k = param[:k]
    initial_state = param[:state].val

    function calc_overlap(state)
        return abs2(state' * initial_state)
    end

    obs = [calc_overlap]

    # diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    # H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]
    # H0 = coupler_H(param) .+ H_dd

    x_op = singlesite_a(d, L, 1) + singlesite_adag(d, L, 1)
    y_op = 1im * (-singlesite_a(d, L, 1) + singlesite_adag(d, L, 1))

    for diss_rate in [0.0, 4]
        diss_op = [sqrt(diss_rate) * singlesite_a(d, L, L)]
        H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]
        H0 = coupler_H(param) .+ H_dd
        res = []
        max_n = 200
        for n in 1:max_n
            function time_evolution()
                H = deepcopy(H0)
                wd = param[:wd].val
                gt = param[:gate_time].val
                p_string = rb_pauli_string(n)
                # println(p_string, " : $n")
                effect! = pauli_string_drive(p_string, wd, gt, x_op, y_op, H, dt; use_rwa = true)
                t = simulation_time(p_string, gt)
                # println("Total time: $t", "Gate time: $gt")
                krylovevolve(initial_state, H, diss_op, dt, t, k, obs...; effect!, save_only_last = true)
            end
            
            r_out = solvetrajectories(time_evolution, traj; paral = :threads)
            push!(res, traj_mean(r_out)[1])
        end
        push!(lines, LineInfo(1:max_n, res, traj, "diss/J = $diss_rate"))
    end

    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "gates", ylabel = "|<psi_0||psi>|^2", copy_code = false, param...)
    addtrajectories(path, lines...)
end


@time f()