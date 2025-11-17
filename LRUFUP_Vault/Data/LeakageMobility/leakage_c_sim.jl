using Distributions
using Statistics


function f()
    width = 7
    height = 7
    n_sites = width * height
    
    frequencies = rand(Normal(4.789996, 0.109104), n_sites)
    anharmonicities = rand(Normal(-0.310393, 0.00548), n_sites)
    J = 0.005 # 5MHz
    Jc = 0.05 # 50MHz

    w1 = 4
    w2 = 4
    wc = wc_idle(w1, w2, J, Jc)
    println("Idle coupler frequency: ", round(wc, digits=4), " GHz")
    J_eff([w1 - 0.3, wc, w2 - 0.3], [Jc, Jc, J])
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

f()