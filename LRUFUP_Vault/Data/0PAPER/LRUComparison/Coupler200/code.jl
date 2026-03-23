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
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")

    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc_below] = ValueInfo(3505 / param[:J12].val, "3.505GHz")

    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    t_10 = 30
    param[:t_10] = ValueInfo(t_10 * 2 * π * param[:J12].val, "$(t_10)μs")

    param[:diss_rate] = 0.08

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
    # n2 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1 + nc, n1, qt_n_all(d, L)]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc_below].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]


    H_below = qt_tunable_coupler_paper(d, w, U, g)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    t_10_list = LinRange(0, param[:t_10].val, param[:t_points])
    
    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    E, psi, Uni = eigenstates(H_below)
    state10 = psi[4]

    state = state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))
    out = mesolve(H_below, state, t_10_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_10_list, real(out.expect[1, :]), 1, "|100>, n1+nc"))

    path = joinpath(@__DIR__, "result_10")
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)

    lines = []
    state = state20
    out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "|200>, n1+nc"))

    path = joinpath(@__DIR__, "result_20")
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)

    lines = []
    state = QuantumObject(tensor(basis(param[:d], 0), basis(param[:d], 1), basis(param[:d], 0)))
    out = mesolve(H_below, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "|010>, n1+nc"))

    path = joinpath(@__DIR__, "result_010")
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
end

@time f()

