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
    N2 = 0
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], 0), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)0$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w1] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:w2] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U1] = ValueInfo(300 / param[:J12].val, "300MHz")
    U2 = 300
    param[:U2] = ValueInfo(U2 / param[:J12].val, "$(U2)MHz")
    param[:Uc] = ValueInfo(300 / param[:J12].val, "300MHz")
    t = 1
    param[:t] = ValueInfo(t * 2 * π * param[:J12].val, "$(t)μs")

    param[:J12] = ValueInfo(1, param[:J12].str_val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------

    d = param[:d]; L = param[:L]
    
    w = [param[:w1].val, param[:wc].val, param[:w2].val]
    U = [param[:U1].val, param[:Uc].val, param[:U2].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]

    H = qt_tunable_coupler_below(w, U, g)
    E, psi, U = eigenstates(H)

    state10 = QuantumObject(tensor(basis(param[:d], 1), basis(param[:d], 0), basis(param[:d], 0)))
    
    for (i, state) in enumerate(psi)
        overlap = abs(state10' * state)^2
        if overlap > 0.9
            println("Eigenvalue $(i): $(E[i])")
            println("Overlap with |100>: ", overlap)
            # println("State: $state")
        end
    end

    println("With leakage")
    state20 = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0), basis(param[:d], 0)))
    n1 = qt_single_site_op(d, L, num(d), 1)
    nc = qt_single_site_op(d, L, num(d), 2)
    n2 = qt_single_site_op(d, L, num(d), 3)
    
    for (i, state) in enumerate(psi)
        overlap = abs(state20' * state)^2
        if overlap > 0.01
            println("--------")
            println("Overlap with |200>: ", overlap)
            println("Eigenvalue $(i): $(E[i])")
            # println("State: $state")
            println("n1: ", expect(n1, state))
            println("nc: ", expect(nc, state))
            println("n2: ", expect(n2, state))
        end
    end
end

f()