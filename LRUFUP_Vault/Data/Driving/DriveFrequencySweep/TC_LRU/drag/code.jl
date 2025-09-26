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
    param[:Uc] = ValueInfo(690 / param[:J12].val, "690MHz")
    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J12].val, "14ns")
    param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:J12].val, "4ns")
    param[:diss_rate] = 3.5

    param[:J12] = ValueInfo(1, "5MHz")
    return param
end

function f()
    param = initial_parameters()
    # drive_frequency_range = LinRange(8.8, 9.4, 20) ./ 5 .+ param[:w].val
    lambda_range = LinRange(0.0, 0.54, 2)
    traj = 100
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    proj_1 = tensor(proj(basis(d, 0)), qeye(d), qeye(d))
    anh_op = qeye(d^L; dims = (3, 3, 3)) - qt_anh_all(d, L) / 2

    c_ops = [sqrt(param[:diss_rate]) * qt_single_site_op(d, L, destroy(d), L)]

    I_op, Q_op = qt_drive_op(d, L; rwa = true)
    
    
    N_MAX = 200
    out = []
    anh_out = []
    for lambda in lambda_range
        pulse, dg_pulse = drag_pulses(param[:gate_time].val, π, lambda, -param[:U].val)
        w = [param[:w].val, param[:wc].val, param[:w].val] .- param[:wd].val
        U = [param[:U].val, -param[:Uc].val, param[:U].val]
        g = [param[:Jc].val, param[:Jc].val, param[:J12].val]
        H = qt_tunable_coupler(w, U, g)
        function time_evolution()
            pauli_string = rb_pauli_string(N_MAX)
            accumulated_phase = 0.0
            state = deepcopy(param[:state].val)
            for p_op in pauli_string
                if p_op == :I
                    Ht = H
                elseif p_op == :X
                    g_part, d_part = qt_drag_gate_to_QobjEvo(accumulated_phase, I_op, Q_op, pulse, dg_pulse)
                    Ht = H + g_part + d_part
                elseif p_op == :Y
                    g_part, d_part = qt_drag_gate_to_QobjEvo(π/2 + accumulated_phase, I_op, Q_op, pulse, dg_pulse)
                    Ht = H + g_part + d_part
                elseif p_op == :Z
                    accumulated_phase += π / 2
                    continue
                end

                t = param[:gate_time].val + param[:gate_repeat_time].val
                t_list = LinRange(0, t, param[:t_points])
                res = mesolve(Ht, state, t_list, c_ops; saveat=[t_list[end]], progress_bar = Val(false))
                state = res.states[1]
            end
            return [real(expect(proj_1, state)), real(expect(anh_op, state))]
        end
        r_out = solvetrajectories(time_evolution, traj; paral = :threads)
        res = traj_mean(r_out)
        push!(out, res[1])
        push!(anh_out, res[2])
    end

    # push!(lines, LineInfo(drive_frequency_range .* 5, out, traj, "proj -> 0"))
    push!(lines, LineInfo(lambda_range, out, traj, "proj -> 0"))
    push!(lines, LineInfo(lambda_range, anh_out, traj, "anh"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "i", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()