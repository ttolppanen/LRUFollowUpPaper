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
    N1 = 0
    Nc = 0
    N2 = 2
    state = QuantumObject(tensor(basis(param[:d], N1), basis(param[:d], Nc), basis(param[:d], N2)))
    param[:state] = ValueInfo(state, "|$(N1)$(Nc)$(N2)>")
    param[:alg] = "toolbox"
    param[:H] = "Tunable Coupler"
    param[:J] = ValueInfo(5, "5MHz")
    param[:Jc] = ValueInfo(50 / param[:J].val, "50MHz")
    param[:w1] = ValueInfo(3700 / param[:J].val, "3.7GHz")
    param[:w2] = ValueInfo(4000 / param[:J].val, "4GHz")
    param[:wc] = ValueInfo(3206 / param[:J].val, "3.206GHz")
    param[:wd] = ValueInfo(3707 / param[:J].val, "3.707GHz")
    param[:U] = ValueInfo(300 / param[:J].val, "300MHz")

    param[:gate_time] = ValueInfo(14 * 10^(-3) * 2 * π * param[:J].val, "14ns")
    param[:gate_repeat_time] = ValueInfo(4 * 10^(-3) * 2 * π * param[:J].val, "4ns")

    param[:J] = ValueInfo(1, "5MHz")
    return param
end

function f()
    param = initial_parameters()
    traj = 200
    # -----------
    # -----------

    lines = []
    d = param[:d]; L = param[:L]

    # proj_1 = tensor(qeye(d), qeye(d), qeye(d)) - qt_single_site_op(d, L, qt_anh_op(d), 1) / 2 # subspace
    # proj_1 = qt_anh_all(d, L) # total anharmonicity
    proj_1 = qt_n_all(d, L) # total excitation number

    label = "Total boson number"

    w = [param[:w1].val, param[:wc].val, param[:w2].val] .- param[:wd].val
    U = [param[:U].val for _ in 1:L]
    g = [param[:Jc].val, param[:Jc].val, param[:J].val]
    H = qt_tunable_coupler(w, U, g)

    I_op, Q_op = qt_drive_op(d, L; rwa = true)
    pulse = gaussian_pulse(param[:gate_time].val, π)
    
    Ns = vcat(2:2:40)
    out = []
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
                res = mesolve(Ht, state, t_list; saveat=[t_list[end]], progress_bar = Val(false))
                state = res.states[1]
            end
            return real(expect(proj_1, state))
        end
        r_out = solvetrajectories(time_evolution, traj; paral = :threads)
        push!(out, traj_mean(r_out))
    end

    push!(lines, LineInfo(Ns, out, traj, label))

    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "i", ylabel = "", copy_code = false, param...)
    # addtrajectories(path, line)
end


@time f()