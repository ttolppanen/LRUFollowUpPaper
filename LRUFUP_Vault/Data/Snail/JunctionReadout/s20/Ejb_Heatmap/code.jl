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


function junction_readout_H(d, E_j, E_c, zpf; D = 15)
    if d > D
        D = d
    end

    Eja = E_j[1]; Ejb = E_j[2]; Ej = E_j[3]
    Eca = E_c[1]; Ecb = E_c[2]; Ec = E_c[3]

    phi_zpf_a = zpf[1]; phi_zpf_b = zpf[2]
    n_zpf_a = 1 / (2 * phi_zpf_a); n_zpf_b = 1 / (2 * phi_zpf_b)
    
    c_op = create(D)
    d_op = destroy(D)

    n_a = -1im * n_zpf_a * (d_op - c_op)
    phi_a = phi_zpf_a * (d_op + c_op)

    Ha = 4 * Eca * n_a^2 - Eja * cos(phi_a)
    λ, psi, U1 = eigenstates(Ha)
    Ua = QuantumObject(U1, dims = D)
    # println(λ[1:3])
    # println((λ[1:3] .- λ[1]))
    # a1 = λ[3] - λ[1]
    
    n_b = -1im * n_zpf_b * (d_op - c_op)
    phi_b = phi_zpf_b * (d_op + c_op)
    Hb = 4 * Ecb * n_b^2 - Ejb / 2 * phi_b^2
    λ, psi, U2 = eigenstates(Hb)
    Ub = QuantumObject(U2, dims = D)
    # println(λ[1:3])
    # println((λ[1:3] .- λ[1]))
    # a2 = λ[2] - λ[1]
    # println("Resonance: $(a1 - a2)")

    U = tensor(Ua, Ub)
    coupling = 8 * Ec * tensor(n_a, n_b) - Ej * (tensor(sin(phi_a), sin(phi_b)) + tensor(cos(phi_a), cos(phi_b)))
    H = U' * (tensor(Ha, qeye(D)) + tensor(qeye(D), Hb) + coupling) * U
    
    H_row = []
    for i in 0:(d-1)
        i_ind = (i * D + 1):(i * D + d)
        push!(H_row, H[(1:d), i_ind])
    end

    H_out = hcat(H_row...)
    for j in 1:(d-1)
        j_ind = (j * D + 1):(j * D + d)
        H_row = []
        for i in 0:(d-1)
            i_ind = (i * D + 1):(i * D + d)
            push!(H_row, H[j_ind, i_ind])
        end
        H_out = vcat(H_out, hcat(H_row...))
    end
    
    H = QuantumObject(H_out, dims = (d, d))
    # λ, psi, U2 = eigenstates(H)
    # E = labeled_eigenvalue_dict(λ, d)
    # println("00 -> 20 = $(E["20"] - E["00"])")
    # # println("00 -> 10 = $(E["10"] - E["00"])")
    # println("00 -> 01 = $(E["01"] - E["00"])")
    return QuantumObject(H_out, dims = (d, d))
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 1000
    param[:d] = 3
    param[:L] = 2

    param[:alg] = "toolbox"
    param[:H] = "j-readout"

    Ej = 1
    param[:Ej] = ValueInfo(1, "$(Ej)GHz")
    param[:Units] = ValueInfo(Ej, "Ej")

    param[:Eja] = ValueInfo(20 / param[:Units].val, "20GHz")
    Ejbs = collect(range(0.1, 1000, 1000)) ./ param[:Units].val
    param[:Ejb] = ValueInfo(Ejbs ./ param[:Units].val, "1-1000GHz")
    param[:Eca] = ValueInfo(0.2 / param[:Units].val, "200MHz")
    param[:Ecb] = ValueInfo(0.045, "45MHz")
    param[:Ec] = ValueInfo(0.0017 / param[:Units].val, "1.7MHz")

    param[:phi_zpf_a] = 0.378
    param[:phi_zpf_b] = 0.165

    t = 10
    param[:t] = ValueInfo(t * 2 * π * param[:Units].val, "$(t)ns")
    # param[:diss_rate] = 0.01
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    t_list = LinRange(0, param[:t].val, param[:t_points])

    n1 = qt_single_site_op(d, L, num(d), 1)
    n2 = qt_single_site_op(d, L, num(d), 2)
    anh1 = qt_single_site_op(d, L, qt_anh_op(d), 1)
    e_ops = [n2]

    E_c = [param[:Eca].val, param[:Ecb].val, param[:Ec].val]
    zpf = [param[:phi_zpf_a], param[:phi_zpf_b]]
    state = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0)))

    populations = zeros(param[:t_points], length(param[:Ejb].val))
    Base.Threads.@threads for i in 1:length(param[:Ejb].val)
        Ejb = param[:Ejb].val[i]

        E_j = [param[:Eja].val, Ejb, param[:Ej].val]
        H = junction_readout_H(d, E_j, E_c, zpf)

        out = mesolve(H, state, t_list, e_ops = e_ops; abstol = 1e-14, reltol = 1e-8, maxiters = 10^7, progress_bar = false)
        populations[:, i] = real.(out.expect[1, :])
    end

    x = param[:Ejb].val .* param[:Units].val
    y = t_list ./ (2 * pi * param[:Units].val)
    heatmap(x, y, populations, xlabel = "Ejb (GHz)", ylabel = "t (ns)", colorbar_title = "<n2>")
    savefig(joinpath(@__DIR__, "heatmap.png"))
end

@time f()