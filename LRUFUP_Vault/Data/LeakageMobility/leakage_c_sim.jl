using Distributions
using Statistics
using Plots
using PlotAndSave
using QuantumTimeEvolution

function f()
    width = 7
    height = 7
    
    J = 5 # 5MHz
    Jc = 50 # 50MHz

    initial_leakage_site = (height ÷ 2 + 1, width ÷ 2 + 1)
    time_step = 0.014 # 14ns
    total_time = 1 # 1mus, 1 round of error correction

    n_traj = 50000
    out = []
    steps = range(0, 1, length=21)
    for step in steps
        frequencies = rand(Normal(4789.996, 109.104 * step), height, width)
        anharmonicities = rand(Normal(-310.393, 5.48 * step), height, width)
        function r_f()
            visited, circuit = simulate_leakage(width, height, frequencies, anharmonicities, J, Jc, initial_leakage_site, time_step, total_time)
            return length(visited)
        end
        lengths = solvetrajectories(r_f, n_traj; paral = :threads)
        push!(out, mean(lengths))
    end

    circuit_size = width * height
    path = joinpath(@__DIR__)
    makeplot(path, LineInfo(steps, out / circuit_size, n_traj, "width = $(width), height = $(height)"); xlabel = "% of Standard Deviation", ylabel = "% of circuit visited", copy_code = false, total_time)
end

function wc_idle(w1, w2, J, Jc) # Maybe missing + J? (+ 5 MHz)?
    wc_plus, wc_minus = solve_wc(w1, w2, -J, Jc) # Minus here to give the idle position under
    return wc_minus
end

function solve_wc(w1, w2, J, Jc)
    return solve_quadratic(1, -(w1 + w2 + Jc^2 / J), w1 * w2 + (w1 + w2) * Jc^2 / (2 * J))
end

function solve_quadratic(a, b, c)
    disc = b^2 - 4*a*c
    if disc < 0
        println("No real solutions. Discriminant = $disc")
        return
    end
    x1 = (-b + sqrt(disc)) / (2*a)
    x2 = (-b - sqrt(disc)) / (2*a)
    # println("x₁ (plus)  = $x1")
    # println("x₂ (minus) = $x2")
    return x1, x2
end

function J_eff(w, g)
    delta_1 = w[1] - w[2]
    delta_2 = w[3] - w[2]
    one_over_delta = 0.5 * (1 / delta_1 + 1 / delta_2)
    j_eff = -g[1] * g[2] * one_over_delta + g[3]
    if abs(j_eff) > (g[1] + g[3])
        return g[1] + g[3]
    end
    return j_eff
end

# 1 is the site we are hopping from, and 2 is the site we are hopping to
function rabi_probability(w1, w2, wc, J, Jc, t)
    Δ = (w1 - w2)
    w = [w1, wc, w2]
    g = [Jc, Jc, J]
    J_e = J_eff(w, g) * 1.62 # magic factor
    JJ = 4 * J_e^2
    omega = sqrt(Δ^2 + JJ)
    return JJ / omega^2 * sin(omega / 2 * t)^2
end

function is_point_inside(point, width, height)
    y = point[1]
    x = point[2]
    return x >= 1 && x <= width && y >= 1 && y <= height
end

function adjacent_sites(site, width, height)
    out = []
    for neighbor in [(-1, 0), (1, 0), (0, -1), (0, 1)]
        new_site = (site[1] + neighbor[1], site[2] + neighbor[2])
        if is_point_inside(new_site, width, height)
            push!(out, new_site)
        end
    end
    return out
end

function random_x_gates!(circuit)
    height, width = size(circuit)
    for y in 1:height
        for x in 1:width
            site_val = circuit[y, x]
            if site_val == 2
                continue
            end
            if rand() < 0.5
                if site_val == 0
                    circuit[y, x] = 1
                else
                    circuit[y, x] = 0
                end
            end
        end
    end
end

function simulate_leakage(width, height, frequencies, anharmonicities, J, Jc, initial_leakage_site, time_step, total_time)
    site = [initial_leakage_site...]
    visited = Set{Tuple{Int, Int}}()
    push!(visited, (site[1], site[2]))
    circuit = zeros(Int, height, width)
    circuit[initial_leakage_site...] = 2

    time_since_jump = 0
    for _ in 1:Int(div(total_time, time_step))
        random_x_gates!(circuit)
        time_since_jump += time_step
        for next_site in adjacent_sites(site, width, height)
            if circuit[next_site...] == 1
                w1 = frequencies[site...]
                w2 = frequencies[next_site...]
                anh1 = anharmonicities[site...]
                anh2 = anharmonicities[next_site...]
                wc = wc_idle(w1, w2, J, Jc)
                prob = rabi_probability(w1 + anh1, w2 + anh2, wc, J, Jc, time_since_jump * 2 * pi)
                if prob > rand()
                    circuit[site...] = 1
                    circuit[next_site...] = 2
                    site = [next_site...]
                    push!(visited, (site[1], site[2]))
                    time_since_jump = 0
                    break
                end
            end
        end
    end
    return visited, circuit
end

@time f()