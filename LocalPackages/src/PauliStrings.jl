using Random

const PAULI_OPS = [:I, :X, :Y, :Z]

export pauli_string
export pauli_string_inverse
export multiply_pauli_gate
export simulation_time
export rb_pauli_string


function pauli_string(n)
    return [rand(PAULI_OPS) for _ in 1:n]
end

function multiply_pauli_gate(a, b)
    if a == :I
        return b, 1
    elseif b == :I
        return a, 1
    elseif a == b
        return :I, 1
    elseif (a == :X && b == :Y)
        return :Z, im
    elseif (a == :Y && b == :X)
        return :Z, -im
    elseif (a == :Y && b == :Z)
        return :X, im
    elseif (a == :Z && b == :Y)
        return :X, -im
    elseif (a == :Z && b == :X)
        return :Y, im
    elseif (a == :X && b == :Z)
        return :Y, -im
    else
        error("Invalid Pauli gates: $a, $b")
    end
end

function pauli_string_inverse(string)
    gate = string[1]
    for i in eachindex(string)
        if i == 1
            continue
        end
        gate, _ = multiply_pauli_gate(string[i], gate)
    end
    return gate # since all pauli gates are self-inverse
end

function simulation_time(pauli_string, gate_time)
    t = 0.0
    for gate in pauli_string
        if gate != :Z
            t += gate_time
        end
    end
    return t
end

function rb_pauli_string(n)
    p_s = pauli_string(n)
    return push!(p_s, pauli_string_inverse(p_s))
end