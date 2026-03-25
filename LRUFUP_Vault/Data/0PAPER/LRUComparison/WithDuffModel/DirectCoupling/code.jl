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
    param[:L] = 2
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "BH"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:U] = ValueInfo(300 / param[:J].val, "300MHz")
    param[:w2] = ValueInfo(param[:w1].val - param[:U].val, "w1 - U")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")
    t_10 = 30
    param[:t_10] = ValueInfo(t_10 * 2 * π * param[:J].val, "$(t_10)μs")
    param[:diss_rate] = 4 # with 0.2 its close to direct coupler one

    param[:J] = ValueInfo(1, param[:J].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]
    t_list = LinRange(0, param[:t].val, param[:t_points])
    t_10_list = LinRange(0, param[:t_10].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L)]
    
    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:w2].val]
    U = [param[:U].val, param[:U].val]
    H = qt_bosehubbard(d, L; w = w, U = U, J = param[:J].val)
    Ht = H

    state = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0)))
    out = mesolve(Ht, state, t_10_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_10_list, real(out.expect[1, :]), 1, "|10>, n1"))

    path = joinpath(@__DIR__, "result_10")
    makeplot(path, lines...; xlabel = "t[μs]", ylabel = "", copy_code = false, param...)

    lines = []
    state = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0)))
    out = mesolve(Ht, state, t_list, c_ops, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "|20>, n1"))

    path = joinpath(@__DIR__, "result_20")
    makeplot(path, lines...; xlabel = "t[μs]", ylabel = "", copy_code = false, param...)
end

# function t_time_from_n(out)
#     x = out.times
#     y = real(out.expect[1, :])
#     pl = plot(x, y)
#     n = length(out.times) # adjust here to make fit better
#     c, slope, f_fit = exp_fit(x[1:n], y[1:n], 1)
#     println(-1 / slope / (2 * π * 5), "μs")
#     plot!(x, f_fit, label = "fit", style = :dash)
# end

@time f()