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
    param[:L] = 2
    N = 1
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)00>")
    param[:alg] = "TB.jl"
    param[:H] = "BH"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:w2] = ValueInfo(3700 / param[:J].val, "2.5GHz")
    param[:U1] = ValueInfo(300 / param[:J].val, "300MHz")
    param[:U2] = ValueInfo(0 / param[:J].val, "0GHz")
    t = 100
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")
    param[:diss_rate] = 0.2
    param[:t1] = 100 * 2 * π * param[:J].val

    param[:J] = ValueInfo(1, param[:J].str_val)
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
    
    c_ops = [sqrt(1 / param[:t1]) * qt_single_site_op(d, L, destroy(d), 1), sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:w2].val]
    U = [param[:U1].val, param[:U2].val]
    H = qt_bosehubbard(d, L; w = w, U = U, J = param[:J].val)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, c_ops, e_ops = e_ops)
    push!(lines, LineInfo(t_list ./ (2 * π * 5), real(out.expect[1, :]), 1, "n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
end


@time f()