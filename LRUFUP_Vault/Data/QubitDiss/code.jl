using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:t] = 500
    param[:d] = 2
    param[:L] = 3
    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = 2), "|100...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(50 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")
    return param
end

function f_original(traj, gammas)
    param = initial_parameters()
    # -----------
    # -----------

    d = param[:d]; L = param[:L]
    dt = param[:dt]; t = param[:t]
    state = param[:state].val

    op = singlesite_n(d, L, 1)
    obs = (state) -> expval(state, op)
    
    lines = []
    for gamma in gammas
        dw = [param[:w].val + param[:dw].val, param[:w].val, param[:w].val - param[:dw].val]
        H = bosehubbard(d, L; w = dw, U = param[:U].val, J = 1)
        diss_op = [sqrt(gamma) * singlesite_a(d, L, L)]
        H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]
        H .+= H_dd
        Uop = exp(-1im * param[:dt] * Matrix(H))

        r_f() = begin
            exactevolve(state, Uop, diss_op, param[:dt], param[:t], obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :none)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "full model, gamma = $gamma"))
    end
    return lines
    # path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end

function f_new_op(traj, gammas)
    param = initial_parameters()
    # -----------

    # -----------

    d = param[:d];
    dt = param[:dt]; t = param[:t]
    state = spzeros(d)
    state[2] = 1.0

    op = nop(d)
    obs = (state) -> expval(state, op)
    
    lines = []
    
    function alpha_squared(dw)
        return (1 + dw^2 - dw * sqrt(2 + dw^2)) / (2 * (2 + dw^2))
    end

    for gamma in gammas
        fixed_gamma = gamma * alpha_squared(param[:dw].val)
        w = param[:w].val + sqrt(2 + param[:dw].val^2)
        H = w * nop(d)
        diss_op = [sqrt(fixed_gamma) * aop(d)]
        H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]
        H .+= H_dd
        Uop = exp(-1im * param[:dt] * Matrix(H))

        r_f() = begin
            exactevolve(state, Uop, diss_op, param[:dt], param[:t], obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :none)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "effective model, gamma = $gamma"))
    end
    return lines
    # path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end

function f_all()
    traj = 3000
    gammas = [0.03, 0.3]

    original_lines = f_original(traj, gammas)
    new_lines = f_new_op(traj, gammas)
    lines = vcat(original_lines, new_lines)
    makeplot(@__DIR__, lines...; xlabel = "t", ylabel = "n1", copy_code = false)
end

@time f_all()