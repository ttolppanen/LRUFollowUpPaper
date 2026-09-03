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
using OrdinaryDiffEq

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 10000
    param[:d] = 3
    param[:L] = 4
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    Jb = 6
    param[:Jb] = ValueInfo(Jb / param[:J12].val, "$(Jb)MHz")
    param[:Jc] = ValueInfo(70 / param[:J12].val, "70MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000/ param[:J12].val, "4GHz")
    param[:wc_off] = ValueInfo(3024.4 / param[:J12].val, "3.0244GHz")
    param[:wc_on] = ValueInfo(3900 / param[:J12].val, "3.9GHz")
    param[:wb] = ValueInfo(3014.45 / param[:J12].val, "3.014.45Hz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    param[:Ub] = ValueInfo(0, "0MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")
    param[:bus_diss_rate] = 24 / param[:J12].val

    # param[:wd] = ValueInfo(4007 / param[:J12].val, "4.007GHz")
    # param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J12].val, "14ns")
    # param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:J12].val, "4ns")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3)
    nb = qt_single_site_op(d, L, num(d), 4)
    e_ops = [n1, nc, n2, nb, qt_n_all(d, L)]
    
    c_ops = [sqrt(param[:bus_diss_rate]) * qt_single_site_op(d, L, destroy(d), 4)]

    w = [param[:w1].val, param[:wc_on].val, param[:w2].val, param[:wb].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val, param[:Ub].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val, param[:Jb].val]


    H = qt_lossy_tunable_coupler(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0), basis(param[:d], 0)))
    state010 = QuantumObject(tensor(basis(param[:d], 0), basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))
    E, psi, Uni = eigenstates(H)
    state10 = psi[5]
    println(round(expect(n1, state10), digits = 3))
    println(round(expect(nc, state10), digits = 3))
    println(round(expect(n2, state10), digits = 3))
    println(round(expect(nb, state10), digits = 3))

    # d_q1, d_q2 = dressed_states(H, n1, n2)

    state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0), basis(param[:d], 0)))

    state = state10
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "s_0 = |100>, <n1>"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "s_0 = |100>, <n2>"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "s_0 = |100>, <n_all>"))

    # state = state010
    # out = mesolve(H, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "s_0 = |010>, <nc>"))
    # push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "s_0 = |010>, <n_all>"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
end

function dressed_states(H, n1, n2)
    E, psi, Uni = eigenstates(H)
    q1 = psi[4]
    q2 = psi[5]
    n = n1 - n2
    h = [q1' * n * q1 q2' * n * q2; q2' * n * q1 q2' * n * q2]
    E, psi, Uni = eigenstates(QuantumObject(h))
    d_q1 = psi[1][1] * q1 + psi[1][2] * q2
    d_q2 = psi[2][1] * q1 + psi[2][2] * q2
    return d_q1, d_q2
end

@time f()

