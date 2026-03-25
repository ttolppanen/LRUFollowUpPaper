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

    Ha = 4 * Eca * n_a^2 - Eja * cos(phi_a) - Ej * cos(phi_a)
    λ, psi, U1 = eigenstates(Ha)
    Ua = QuantumObject(U1, dims = D)
    # println(λ[1:3])
    println("transmon")
    println((λ[1:3] .- λ[1]))
    a1 = λ[3] - λ[1]
    
    n_b = -1im * n_zpf_b * (d_op - c_op)
    phi_b = phi_zpf_b * (d_op + c_op)
    Hb = 4 * Ecb * n_b^2 + Ejb / 2 * phi_b^2 - Ej * cos(phi_b)
    λ, psi, U2 = eigenstates(Hb)
    Ub = QuantumObject(U2, dims = D)
    # println(λ[1:3])
    println("resonator")
    println((λ[1:3] .- λ[1]))
    normalized_E = λ[1:3] .- λ[1]
    println(normalized_E[3] / normalized_E[2])
    a2 = λ[3] - λ[1]
    println("Resonance: $(a1 - a2)")

    U = tensor(Ua, Ub)
    coupling = 8 * Ec * tensor(n_a, n_b) - Ej * (tensor(sin(phi_a), sin(phi_b)) + tensor(cos(phi_a), cos(phi_b)) - tensor(cos(phi_a), qeye(D)) - tensor(qeye(D), cos(phi_b)))
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
    param[:Ejb] = ValueInfo(100.46 / param[:Units].val, "100GHz")
    param[:Eca] = ValueInfo(0.2 / param[:Units].val, "200MHz")
    param[:Ecb] = ValueInfo(0.037 / param[:Units].val, "37MHz")
    param[:Ec] = ValueInfo(0.00178 / param[:Units].val, "1.7MHz")

    param[:phi_zpf_a] = 0.378
    param[:phi_zpf_b] = 0.165

    t = 300
    param[:t] = ValueInfo(t * 2 * π * param[:Units].val, "$(t)ns")
    param[:diss_rate] = 0.02 / param[:Units].val
    param[:t1] = 1 / (50 * 1e3 * 2 * pi * param[:Units].val)
    return param
end

function f()
    param = initial_parameters()
    # -----------
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    

    n1 = qt_single_site_op(d, L, num(d), 1)
    n2 = qt_single_site_op(d, L, num(d), 2)
    anh1 = qt_single_site_op(d, L, qt_anh_op(d), 1)
    anh2 = qt_single_site_op(d, L, qt_anh_op(d), 2)
    e_ops = [n1, n2, anh1]
    
    c_ops_both = [sqrt(param[:t1]) * qt_single_site_op(d, L, destroy(d), 1), sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]
    c_ops_T1 = [deepcopy(c_ops_both[1])]

    E_j = [param[:Eja].val, param[:Ejb].val, param[:Ej].val]
    E_c = [param[:Eca].val, param[:Ecb].val, param[:Ec].val]
    zpf = [param[:phi_zpf_a], param[:phi_zpf_b]]
    H = junction_readout_H(d, E_j, E_c, zpf)

    t_list = LinRange(0, param[:t].val, param[:t_points])

    state = QuantumObject(tensor(basis(param[:d], 2), basis(param[:d], 0)))
    # out = mesolve(H, state, t_list, c_ops_T1, e_ops = e_ops; alg = Vern9(), abstol = 1e-14, reltol = 1e-8, maxiters = 10^7)
    # push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "T1 only, n1"))
    # push!(lines, LineInfo(t_list, real(out.expect[2, :]), 1, "|20>, n2"))
    # push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "|20>, anh1"))

    out = mesolve(H, state, t_list, c_ops_both, e_ops = e_ops; alg = Vern9(), abstol = 1e-14, reltol = 1e-8, maxiters = 10^7)
    push!(lines, LineInfo(t_list, real(out.expect[1, :]), 1, "T1 and LRU dissipation, n1"))
    push!(lines, LineInfo(t_list, real(out.expect[3, :]), 1, "|20>, anh1"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)

    # drive_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    # drive_op = drive_op / 2
    # pulse = gaussian_pulse(param[:gate_time].val, π)
    # Ht = H + QobjEvo(drive_op, (p, t) -> pulse(t))
    # pl = plot(lines[1].x, lines[1].y, ylim = (0, 1.1))
    # plot!(pl, lines[1].x, t -> exp(-t / (2 * π * 20.1) / 30))
    # t_time_from_n(t_list, real(out.expect[1, :]), 10 * 1e9)
end

@time f()