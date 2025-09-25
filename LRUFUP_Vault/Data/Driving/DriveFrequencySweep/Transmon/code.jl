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
    param[:L] = 1
    N = 0
    state = QuantumObject(basis(param[:d], N))
    param[:state] = ValueInfo(state, "|$(N)>")
    param[:alg] = "toolbox"
    param[:H] = "BH"

    param[:U] = ValueInfo(300, "300MHz")
    param[:w] = ValueInfo(4000 / param[:U].val, "4GHz")
    param[:wd] = ValueInfo(4004 / param[:U].val, "4004GHz")

    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:U].val, "20ns")
    param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:U].val, "4ns")

    param[:U] = ValueInfo(1, param[:U].str_val)
    return param
end

function f()
    param = initial_parameters()
    drive_frequency_range = LinRange(-20, 20, 10) ./ 300 .+ param[:w].val  # Example: adjust as needed
    traj = 500

    lines = []
    d = param[:d]; L = param[:L]

    proj_1 = qeye(d) - qt_anh_op(d) / 2

    I_op, Q_op = qt_drive_op(d, L; rwa = true)
    pulse = gaussian_pulse(param[:gate_time].val, π)

    N_MAX = 10  # Fixed Pauli string length for each frequency
    out = []
    for wd in drive_frequency_range
        w = param[:w].val - wd
        H = qt_bosehubbard(d, L; w, U = param[:U].val, J = 0.0)

        function time_evolution()
            pauli_string = rb_pauli_string(N_MAX)
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
                res = mesolve(Ht, state, t_list; saveat=[t_list[end]], progress_bar = Val(false))
                state = res.states[1]
            end
            return real(expect(proj_1, state))
        end
        r_out = solvetrajectories(time_evolution, traj; paral = :threads)
        push!(out, traj_mean(r_out))
    end

    push!(lines, LineInfo(drive_frequency_range .* 300, out, traj, "|2> population"))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "Drive frequency (MHz)", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()