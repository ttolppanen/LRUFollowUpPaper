using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox

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
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    state = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|100>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wd] = ValueInfo(4009.1 / param[:J12].val, "4.02GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U] = ValueInfo(250 / param[:J12].val, "250MHz")
    param[:Uc] = ValueInfo(732 / param[:J12].val, "740MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:gate_time] = ValueInfo(20 * 10^(-3) * 2 * π * param[:J12].val, "20ns")
    param[:gate_area] = ValueInfo(π, "π")
    param[:gate_repeat_time] = ValueInfo(0.01 * 2 * π * param[:J12].val, "1μs")
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
    state = param[:state].val

    e_ops = [tensor(num(d), qeye(d), qeye(d))]
    
    # diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    # H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]

    H = coupler_H(param)
    H = QuantumObject(H, dims = (d, d, d))
    # H = coupler_H(param) .+ H_dd
    
    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()