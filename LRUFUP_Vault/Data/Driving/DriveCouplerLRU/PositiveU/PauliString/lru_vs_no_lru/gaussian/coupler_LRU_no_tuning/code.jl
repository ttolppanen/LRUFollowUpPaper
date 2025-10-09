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

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:t_points] = 2
    param[:d] = 3
    param[:L] = 3
    N = 0
    state = QuantumObject(tensor(basis(param[:d], N), basis(param[:d], 0), basis(param[:d], 0)))
    param[:state] = ValueInfo(state, "|$(N)00>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J12] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J12].val, "50MHz")
    param[:w] = ValueInfo(4000 / param[:J12].val, "4GHz")
    param[:wd] = ValueInfo(4009.1 / param[:J12].val, "4.0091GHz")
    param[:wc] = ValueInfo(3505 / param[:J12].val, "3.505GHz")
    param[:U] = ValueInfo(300 / param[:J12].val, "300MHz")
    param[:Uc] = ValueInfo(param[:U].val, "U")
    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J12].val, "14ns")
    param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:J12].val, "4ns")
    param[:diss_rate] = 0.08

    param[:J12] = ValueInfo(1, "5MHz")
    return param
end

function f()
    param = initial_parameters()
    traj = 100
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    anh = tensor(qeye(d), qeye(d), qeye(d)) - qt_anh_all(d, L) / 2
    proj_1 = tensor(proj(basis(d, 0)), qeye(d), qeye(d))

    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    w = [param[:w].val, param[:wc].val, param[:w].val] .- param[:wd].val
    U = [param[:U].val, param[:Uc].val, param[:U].val]
    g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
    H = qt_tunable_coupler(w, U, g)

    I_op, Q_op = qt_drive_op(d, L; rwa = true)
    pulse = gaussian_pulse(param[:gate_time].val, π)
    
    Ns = vcat(1, 10:10:500)
    out = []
    anh_out = []
    for n in Ns
        function time_evolution()
            pauli_string = rb_pauli_string(n)
            accumulated_phase = 0.0
            state = deepcopy(param[:state].val)
            for p_op in pauli_string
                if p_op == :I
                    Ht = H
                elseif p_op == :X
                    Ht = H + qt_x_gate_to_QobjEvo(accumulated_phase, I_op, Q_op, pulse)
                elseif p_op == :Y
                    Ht = H + qt_y_gate_to_QobjEvo(accumulated_phase, I_op, Q_op, pulse)
                elseif p_op == :Z
                    accumulated_phase += π / 2
                    continue
                end

                t = param[:gate_time].val + param[:gate_repeat_time].val
                t_list = LinRange(0, t, param[:t_points])
                res = mesolve(Ht, state, t_list, c_ops; saveat=[t_list[end]], progress_bar = Val(false))
                state = res.states[1]
            end
            return [real(expect(proj_1, state)), real(expect(anh, state))]
        end
        r_out = solvetrajectories(time_evolution, traj; paral = :threads)
        r_out = traj_mean(r_out)
        push!(out, r_out[1])
        push!(anh_out, r_out[2])
    end

    push!(lines, LineInfo(Ns, out, traj, "<Ψ|0><0|Ψ>"))
    push!(lines, LineInfo(Ns, anh_out, traj, "Qubit subspace population"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "i", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()