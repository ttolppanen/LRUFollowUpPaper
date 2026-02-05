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
    N1 = 1
    N2 = 0
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)0$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    w2s = vcat(collect(range(3995, 4000, length=50)), collect(range(4000, 4005, length=50))) ./ param[:J12].val
    param[:w2] = ValueInfo(w2s, "(3.995-4.005)GHz")
    wcs = collect(range(3000, 4500, length=1200)) ./ param[:J12].val
    param[:wc] = ValueInfo(wcs, "(3-4.5)GHz")
    # param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    # param[:wc_below] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    # param[:wc_above] = ValueInfo(4495 / param[:J12].val, "4.495GHz")
    U = 300
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t_below = 1
    param[:t_below] = ValueInfo(t_below * 2 * π * param[:J12].val, "$(t_below)μs")
    t_above = 15
    param[:t_above] = ValueInfo(t_above * 2 * π * param[:J12].val, "$(t_above)μs")
    # param[:diss_rate] = 0.0

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
    t_list_below = LinRange(0, param[:t_below].val, param[:t_points])
    t_list_above = LinRange(0, param[:t_above].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3) / 2
    e_ops = [n2]
    
    # c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))

    max_amplitude_below = []
    for w2 in param[:w2].val
        wc = w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = false)
        w = [param[:w1].val, wc, w2]
        U = [param[:U1].val, param[:Uc].val, param[:U2].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

        H = qt_tunable_coupler_paper(w, U, g)

        state = state20
        out = mesolve(H, state, t_list_below, e_ops = e_ops, progress_bar= false)
        push!(max_amplitude_below, maximum(real(out.expect[1, :])))
    end
    push!(lines, LineInfo(param[:w2].val .* 5, max_amplitude_below, 1, "g = 0, max(<n2 / 2>)"))

    # max_amplitude_above = []
    # for w2 in param[:w2].val
    #     println(w2)
    #     wc = w_off_position(param[:w1].val, w2, param[:Jc].val, param[:Jc].val, param[:J12].val; is_above = true)
    #     w = [param[:w1].val, wc, w2]
    #     U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    #     g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    #     H = qt_tunable_coupler_paper_positive(w, U, g)

    #     state = state20
    #     out = mesolve(H, state, t_list_above, e_ops = e_ops, progress_bar= false, maxiters = 10^8)
    #     push!(max_amplitude_above, maximum(real(out.expect[3, :])))
    # end
    # push!(lines, LineInfo(param[:w2].val .* 5, max_amplitude_above, 1, "above, max(<n2>)"))

    max_amplitude_zz_off = []
    for w2 in param[:w2].val
        U = [param[:U1].val, param[:Uc].val, param[:U2].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
        wc = zz_zero_point(param[:w1].val, w2, param[:wc].val, U, g)
        w = [param[:w1].val, maximum(wc), w2]

        H = qt_tunable_coupler_paper(w, U, g)

        state = state20
        out = mesolve(H, state, t_list_below, e_ops = e_ops, progress_bar= false)
        push!(max_amplitude_zz_off, maximum(real(out.expect[1, :])))
    end
    push!(lines, LineInfo(param[:w2].val .* 5, max_amplitude_zz_off, 1, "zz = 0, max(<n2 / 2>)"))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "MHz", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()