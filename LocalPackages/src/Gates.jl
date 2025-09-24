# using Integrals

export gaussian_pulse
export drag_pulses
export drive_hamiltonian
export rand_drive_hamiltonian
export pauli_string_drive
export qt_x_gate_to_QobjEvo
export qt_y_gate_to_QobjEvo
export qt_drag_gate_to_QobjEvo
export qt_drive_op


mutable struct PauliPulseInfo
    time::Float64
    phase::Float64
    pauli_string::Vector{Symbol}
    i::Int # index of the pauli string
end


function gaussian_pulse(gate_time, area)
    f(t, p) = exp(-(t - gate_time / 2)^2 / (2 * gate_time^2)) - exp(-(gate_time / 2)^2 / (2 * gate_time^2))
    prob = IntegralProblem(f, (0.0, gate_time))
    sol = solve(prob, QuadGKJL())
    pulse_area = area / sol.u

    function pulse(t)
        if t < 0 || t > gate_time
            return 0.0
        else
            return pulse_area * f(t, 0)
        end
    end
    return pulse
end

function drag_pulses(gate_time, area, lambda, delta)
    g_pulse = gaussian_pulse(gate_time, area)
    dg_pulse(t) = -(lambda / delta) * (t - gate_time / 2) / (gate_time^2) * g_pulse(t)
    return g_pulse, dg_pulse
end

function drive_hamiltonian(w_d, gate_time, area, drive_op, H, dt; repeat_time = -1, use_rwa = false)
    total_time = [0.0]
    H0 = deepcopy(H)

    pulse = gaussian_pulse(gate_time, area)
    
    function Ht(state)
        total_time[1] += dt
        if repeat_time >= 0.0
            t = total_time[1] % (gate_time + repeat_time) 
        else
            t = total_time[1]
        end

        H .= H0
        pulse_strength = pulse(t)
        if pulse_strength > 0.0
            if use_rwa
                H .+= -pulse_strength / 2 * drive_op
            else
                H .+= -pulse_strength * cos(w_d * t) * drive_op
            end
        end
        return state
    end
    return Ht
end

function rand_drive_hamiltonian(w_d, gate_time, area, I_op, Q_op, H, dt; repeat_time = -1, use_rwa = false)
    total_time = [0.0]
    H0 = deepcopy(H)
    drive_op = deepcopy(I_op) # just reserving memory

    pulse = gaussian_pulse(gate_time, area)
    
    function Ht(state)
        total_time[1] += dt
        if repeat_time >= 0.0
            if total_time[1] >= gate_time + repeat_time
                rand_phase = 2 * π * rand()
                drive_op .= cos(rand_phase) * I_op + sin(rand_phase) * Q_op
                total_time[1] -= gate_time + repeat_time
            end
        end
        t = total_time[1]

        H .= H0
        pulse_strength = pulse(t)
        if pulse_strength > 0.0
            if use_rwa
                H .+= -pulse_strength / 2 * drive_op
            else
                H .+= -pulse_strength * cos(w_d * t) * drive_op
            end
        end
        return state
    end
    return Ht
end

# area is assumed to be π, and also the inverse of the string is added to the end
function pauli_string_drive(p_string, w_d, gate_time, I_op, Q_op, H, dt; use_rwa = false)
    H0 = deepcopy(H)
    drive_op = deepcopy(I_op) # just reserving memory
    pulse_info = PauliPulseInfo(0.0, 0.0, p_string, 1)
    drive_op_from_pauli_gate!(drive_op, I_op, Q_op, pulse_info)
    pulse_envelope = gaussian_pulse(gate_time, π)

    function Ht(state)
        pulse_info.time += dt
        if pulse_info.time >= gate_time
            drive_op_from_pauli_gate!(drive_op, I_op, Q_op, pulse_info)
            pulse_info.time -= gate_time
        end
        t = pulse_info.time
        H .= H0
        
        if pulse_info.pauli_string[pulse_info.i - 1] == :I
            return state
        end

        pulse_strength = pulse_envelope(t)
        if pulse_strength > 0.0
            if use_rwa
                H .+= -pulse_strength / 2 * drive_op
            else
                H .+= -pulse_strength * cos(w_d * t) * drive_op
            end
        end
        return state
    end
    return Ht
end

function drive_op_from_pauli_gate!(drive_op, I_op, Q_op, pulse_info)
    if pulse_info.i > length(pulse_info.pauli_string)
        return
    end
    gate = pulse_info.pauli_string[pulse_info.i]
    pulse_info.i += 1
    if gate == :I
    elseif gate == :X
        drive_op .= cos(0.0 + pulse_info.phase) * I_op + sin(0.0 + pulse_info.phase) * Q_op
    elseif gate == :Y
        drive_op .= cos(π/2 + pulse_info.phase) * I_op + sin(π/2 + pulse_info.phase) * Q_op
    elseif gate == :Z
        pulse_info.phase = mod(pulse_info.phase + π, 2π)
        drive_op_from_pauli_gate!(drive_op, I_op, Q_op, pulse_info)
    else
        error("Invalid Pauli gate: $gate")
    end
end

#     elseif gate == :Z
#         accumulated_phase = mod(pulse_info.phase + π / 2, 2π)
# test for dynamic dispatch stuff
function qt_x_gate_to_QobjEvo(accumulated_phase, I_op::QuantumObject, Q_op::QuantumObject, pi_pulse)
    return QobjEvo(cos(0.0 + accumulated_phase) * I_op + sin(0.0 + accumulated_phase) * Q_op, (p, t) -> pi_pulse(t))
end

function qt_y_gate_to_QobjEvo(accumulated_phase, I_op::QuantumObject, Q_op::QuantumObject, pi_pulse)
    return QobjEvo(cos(π/2 + accumulated_phase) * I_op + sin(π/2 + accumulated_phase) * Q_op, (p, t) -> pi_pulse(t))
end

function qt_drag_gate_to_QobjEvo(phase, I_op::QuantumObject, Q_op::QuantumObject, pi_pulse, pi_dg_pulse)
    gaussian_drive = QobjEvo(cos(phase) * I_op + sin(phase) * Q_op, (p, t) -> pi_pulse(t))
    drag_drive = QobjEvo(-sin(phase) * I_op + cos(phase) * Q_op, (p, t) -> pi_dg_pulse(t)) # minus here might not be right?
    return gaussian_drive, drag_drive
end

function qt_drive_op(d, L; rwa = false)
    I_op = qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1)
    Q_op = im * (-qt_single_site_op(d, L, destroy(d), 1) + qt_single_site_op(d, L, create(d), 1))
    if rwa
        I_op = I_op / 2
        Q_op = Q_op / 2
    end
    return I_op, Q_op
end