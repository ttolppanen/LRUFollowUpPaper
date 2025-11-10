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
using QuantumTimes

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    N = 2
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)00>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wd] = ValueInfo(4009.1 / param[:J12].val, "4.0091GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    param[:Uc] = ValueInfo(-600 / param[:J12].val, "300MHz")
    t = 0.1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J12].val, "14ns")
    param[:gate_area] = ValueInfo(π, "π")
    param[:gate_repeat_time] = ValueInfo(0.01 * 2 * π * param[:J12].val, "1μs")
    param[:diss_rate] = 0.00

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

    n1 = qt_single_site_op(d, L, num(d), 1)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L)]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w].val, param[:wc].val, param[:w].val] .- param[:wd].val
    U = [param[:U].val, param[:Uc].val, param[:U].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
    H = qt_tunable_coupler(w, U, g)
    drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    drive_op = drive_op / 2
    pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))
    Ht = H

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(Ht, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

function t_time_from_n(times, population)
    x = out.times
    y = real(out.expect[1, :])
    pl = plot(x, y)
    n = length(out.times) # adjust here to make fit better
    c, slope, f_fit = exp_fit(x[1:n], y[1:n], 2)
    println(-1 / slope / (2 * π * 5 * 1e6 * 1e-9), "ns")
    plot!(x, f_fit)
end

@time f()