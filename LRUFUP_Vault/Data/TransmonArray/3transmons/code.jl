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
    N1 = 2
    N2 = 1
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], N2), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w1] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:w2] = ValueInfo(3960 / param[:J].val, "3.960GHz")
    param[:w3] = ValueInfo(3995 / param[:J].val, "3.995GHz")
    U = 300
    param[:U] = ValueInfo(U / param[:J].val, "$(U)MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J].val, "$(t)μs")
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
    n2 = qt_single_site_op(d, L, num(d), 2)
    n3 = qt_single_site_op(d, L, num(d), 3)
    e_ops = [n1, qt_anh_all(d, L), qt_n_all(d, L), n2, n3]

    w = [param[:w1].val, param[:w2].val, param[:w3].val]
    H = qt_bosehubbard(d, L; J = param[:J].val, U = param[:U].val, w)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    out = mesolve(H, state, t_list, e_ops = e_ops)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "n1"))

    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "n_all"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "n2"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "n3"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "anh_all"))


    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()