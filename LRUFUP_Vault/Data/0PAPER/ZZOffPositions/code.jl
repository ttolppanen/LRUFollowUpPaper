using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using OrderedCollections
using SparseArrays
using LocalPackages
using QuantumToolbox


function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 3

    param[:J12] = ValueInfo(3.7, "5MHz")
    param[:Jc] = ValueInfo(52 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4100 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(3890 / param[:J12].val, "4GHz")
    wcs = collect(range(3100, 3200, length=100)) ./ param[:J12].val
    param[:wc] = ValueInfo(wcs, "(3-4)GHz")
    U = 216
    param[:U1] = ValueInfo(U / param[:J12].val, "$(U)MHz")
    U2 = U
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    Uc = U
    param[:Uc] = ValueInfo(Uc / param[:J12].val, "$(Uc)MHz")
    t = 10

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    zzs = []
    comparison_states = []

    println("Calculating zz coupling...")
    println(3.7 * zz_zero_point(param[:w1].val, param[:w2].val, param[:wc].val,
        [param[:U1].val, param[:Uc].val, param[:U2].val],
        [param[:Jc].val, param[:Jc].val, param[:J12].val]))
    for wc in param[:wc].val
        if isempty(comparison_states)
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 1))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 1))))
        end
        w = [param[:w1].val, wc, param[:w2].val]
        U = [param[:U1].val, param[:Uc].val, param[:U2].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
        H = qt_tunable_coupler_paper(w, U, g)

        E, psi, Uni = eigenstates(H)
        dressed_states_i = logical_state_index(comparison_states, psi, d)
        push!(zzs, zz_coupling_strength(dressed_states_i, E, d))
        comparison_states = psi[collect(dressed_states_i)]
    end
    push!(lines, LineInfo(param[:wc].val .* 3.7, zzs, 1, ""))
    
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "zz coupling", copy_code = false, param...)
end

function zz_zero_point(w1, w2, wcs, Us, gs)
    zzs = []
    comparison_states = []
    d = 3
    for wc in wcs
        if isempty(comparison_states)
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 0))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 0), basis(d, 0), basis(d, 1))))
            push!(comparison_states, QuantumObject(tensor(basis(d, 1), basis(d, 0), basis(d, 1))))
        end
        w = [w1, wc, w2]
        U = Us
        g = gs
        H = qt_tunable_coupler_paper(w, U, g)

        E, psi, Uni = eigenstates(H)
        dressed_states_i = logical_state_index(comparison_states, psi, d)
        push!(zzs, zz_coupling_strength(dressed_states_i, E, d))
        comparison_states = psi[collect(dressed_states_i)]
    end
    wc_zeros = []
    prev_zz_sign = sign(zzs[1])
    for (i, zz) in enumerate(zzs)
        if sign(zz) != prev_zz_sign
            push!(wc_zeros, wcs[i])
            prev_zz_sign = sign(zz)
        end
    end
    return maximum(wc_zeros)
end

f()