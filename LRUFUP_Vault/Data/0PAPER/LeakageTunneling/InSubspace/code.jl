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
    param[:d] = 5
    param[:L] = 5
    param[:alg] = "toolbox Vern9"
    param[:abstol] = 1e-14
    param[:reltol] = 1e-8
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    w = 3719.50
    param[:w1] = ValueInfo(w / param[:J12].val, "3.75GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w3] = ValueInfo(w / param[:J12].val, "3.75GHz")
    wc = zz_0_with_w_4GHZ(w)
    param[:wc1] = ValueInfo(wc / param[:J12].val, "3.162762GHz")
    param[:wc2] = ValueInfo(wc / param[:J12].val, "3.162762GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    # U = [param[:U].val, param[:U].val, param[:U].val]
    # gs = [param[:Jc].val, param[:Jc].val, 1]
    # wcs = collect(range(2000, 4200, length=1000)) ./ param[:J12].val
    # wcs = zz_zero_point(param[:w1].val, param[:w2].val, wcs, U, gs)
    # display(wcs .* param[:J12].val)

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    indices = total_boson_number_subspace_indices(d, L)
    project_to_subspace(m) = subspace_split(m, indices)

    n1 = singlesite_n(d, L, 1); n1 = project_to_subspace(n1)
    n2 = singlesite_n(d, L, 3); n2 = project_to_subspace(n2)
    n3 = singlesite_n(d, L, 5); n3 = project_to_subspace(n3)
    nc1 = singlesite_n(d, L, 2); nc1 = project_to_subspace(nc1)
    nc2 = singlesite_n(d, L, 4); nc2 = project_to_subspace(nc2)

    # e_ops = [n1, n2, n3, nc1, nc2, qt_n_all(d, L)]
    
    # # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w1].val, param[:wc1].val, param[:w2].val, param[:wc2].val, param[:w3].val]
    U = [param[:U].val, param[:U].val, param[:U].val, param[:U].val, param[:U].val]
    g = [param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:Jc].val, param[:J12].val]

    H = tunnel_tunable_coupler_paper(d, w, U, g); H = project_to_subspace(H)

    t_list = LinRange(0, param[:t].val, param[:t_points])
    
    state200 = productstate(d, [2, 0, 0, 0, 0]); initial_id_200 = find_subspace(state200, indices); state200 = project_to_subspace(state200)
    state201 = productstate(d, [2, 0, 0, 0, 1]); initial_id_201 = find_subspace(state201, indices); state201 = project_to_subspace(state201)
    state210 = productstate(d, [2, 0, 1, 0, 0]); initial_id_210 = find_subspace(state210, indices); state210 = project_to_subspace(state210)
    state211 = productstate(d, [2, 0, 1, 0, 1]); initial_id_211 = find_subspace(state211, indices); state211 = project_to_subspace(state211)
 
    # E, psi, Uni = eigenstates(H)
    # state10 = psi[4]

    # CHANGE HERE
    id = initial_id_200
    state = state200
    label_s = "200"
    
    state = QuantumObject(state[id])
    H = QuantumObject(H[id])
    e_ops = [n1, n2, n3, nc1, nc2]
    e_ops = [QuantumObject(e[id]) for e in e_ops]
    out = mesolve(H, state, t_list, e_ops = e_ops; alg = Vern9(), abstol = param[:abstol], reltol = param[:reltol], maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "s_0 = |$label_s>, <n1>"))
    push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "s_0 = |$label_s>, <n2>"))
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "s_0 = |$label_s>, <n3>"))
    push!(lines, LineInfo(t_list, real(out.expect[4, :]), 1, "s_0 = |$label_s>, <nc1>"))
    push!(lines, LineInfo(t_list, real(out.expect[5, :]), 1, "s_0 = |$label_s>, <nc2>"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
end

# initial_parameters()

@time f()