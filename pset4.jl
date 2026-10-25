
# BEEM138 Macroeconomic Theory I — Problem Set 4
# Neoclassical growth model: sequential shooting algorithm

using Printf
using Plots

# =========================================================
# 1. PARAMETERS GIVEN IN THE ASSIGNMENT
# =========================================================
A = 1.0
alpha = 0.3
sigma = 2.0
rho = 0.05
delta = 0.05
beta = 1 / (1 + rho)
T = 60

f(k) = A * k^alpha
fp(k) = A * alpha * k^(alpha - 1)

# =========================================================
# 2. PART (a): ANALYTICAL STEADY STATE
# =========================================================
# At the steady state:
# 1 = beta * (f'(k*) + 1 - delta)
# Hence f'(k*) = rho + delta.

kstar = (A * alpha / (rho + delta))^(1 / (1 - alpha))
cstar = f(kstar) - delta * kstar
k0 = 0.5 * kstar

# =========================================================
# 3. FORWARD SIMULATION FOR A GIVEN INITIAL CONSUMPTION
# =========================================================
function simulate(c0, k0, T)
    capital = fill(NaN, T + 1)
    consumption = fill(NaN, T + 1)

    capital[1] = k0
    consumption[1] = c0

    for t in 1:T
        k = capital[t]
        c = consumption[t]

        # Resource constraint
        knew = f(k) + (1 - delta) * k - c

        # Reject paths with non-positive or invalid capital
        if !isfinite(knew) || knew <= 0
            return nothing
        end

        # Euler equation
        euler_factor = beta * (fp(knew) + 1 - delta)

        if !isfinite(euler_factor) || euler_factor <= 0
            return nothing
        end

        cnew = c * euler_factor^(1 / sigma)

        if !isfinite(cnew) || cnew <= 0
            return nothing
        end

        capital[t + 1] = knew
        consumption[t + 1] = cnew
    end

    return (k = capital, c = consumption)
end

# =========================================================
# 4. PART (b): SHOOTING ALGORITHM
# =========================================================
# 4. PART (b): ROBUST SHOOTING ALGORITHM

function shooting_gap(c0, k0, T)
    kt = k0
    ct = c0

    for t in 1:T
        knew = f(kt) + (1 - delta) * kt - ct

        if !isfinite(knew) || knew > 1e12
            return Inf
        end

        if knew <= 0
            return -Inf
        end

        factor = beta * (fp(knew) + 1 - delta)

        if !isfinite(factor) || factor <= 0
            return -Inf
        end

        cnew = ct * factor^(1 / sigma)

        if !isfinite(cnew) || cnew <= 0
            return -Inf
        end

        kt = knew
        ct = cnew
    end

    return kt - kstar
end

function find_initial_consumption(k0, T)
    low = 1e-8
    high = f(k0) + (1 - delta) * k0 - 1e-8

    gap_low = shooting_gap(low, k0, T)
    gap_high = shooting_gap(high, k0, T)

    if !(gap_low > 0 && gap_high < 0)
        error(
            "Shooting endpoints do not bracket the solution. " *
            "gap_low = $gap_low, gap_high = $gap_high. " *
            "Check the model equations and parameters."
        )
    end

    for iteration in 1:200
        mid = (low + high) / 2
        gap_mid = shooting_gap(mid, k0, T)

        if gap_mid > 0
            low = mid
        else
            high = mid
        end

        if abs(high - low) < 1e-13
            break
        end
    end

    return (low + high) / 2
end

c0 = find_initial_consumption(k0, T)
path = simulate(c0, k0, T)

if path === nothing
    error("The final simulated path is infeasible.")
end

k = path.k
c = path.c
time = 0:T

println()
println("PART (b): SHOOTING ALGORITHM")
@printf("Initial capital k0 = %.5f\n", k0)
@printf("Initial consumption c0 = %.5f\n", c0)
@printf("Terminal capital kT = %.5f\n", k[end])
@printf("Steady-state capital k* = %.5f\n", kstar)

# =========================================================
# 5. PART (c): HALF-LIFE
# =========================================================
halfway_k = k0 + 0.5 * (kstar - k0)
half_life = findfirst(x -> x >= halfway_k, k)

# Julia array index 1 represents time t = 0.
half_life_period = isnothing(half_life) ? nothing : half_life - 1

# =========================================================
# 6. SAVE NUMERICAL RESULTS
# =========================================================
mkpath("results")
mkpath("figures")

open("results/summary.txt", "w") do io
    println(io, "BEEM138 Problem Set 4 — Results")
    println(io, "PART (a): Steady state")
    @printf(io, "k* = %.10f\n", kstar)
    @printf(io, "c* = %.10f\n\n", cstar)

    println(io, "PART (b): Shooting algorithm")
    @printf(io, "k0 = %.10f\n", k0)
    @printf(io, "c0 = %.12f\n", c0)
    @printf(io, "k_T = %.10f\n", k[end])
    @printf(io, "c_T = %.10f\n", c[end])
    @printf(io, "Terminal capital gap = %.10f\n", k[end] - kstar)
    @printf(io, "Terminal consumption gap = %.10f\n\n", c[end] - cstar)

    println(io, "PART (c): Half-life")
    @printf(io, "Halfway capital target = %.10f\n", halfway_k)
    if isnothing(half_life_period)
        println(io, "Halfway point not reached within T periods.")
    else
        @printf(io, "Discrete half-life = %d periods\n", half_life_period)
    end

    println(io, "\nPART (d): Policy path saved to results/policy_function.csv")
end

open("results/time_series.csv", "w") do io
    println(io, "t,k_t,c_t")
    for i in eachindex(time)
        @printf(io, "%d,%.12f,%.12f\n", time[i], k[i], c[i])
    end
end

open("results/policy_function.csv", "w") do io
    println(io, "k_t,k_t_plus_1")
    for i in 1:T
        @printf(io, "%.12f,%.12f\n", k[i], k[i + 1])
    end
end

# =========================================================
# 7. PART (b): SADDLE PATH AND DELTA K = 0 CURVE
# =========================================================
kgrid = range(0.1, 1.35 * kstar, length=500)
c_delta_k_zero = [f(x) - delta * x for x in kgrid]

p1 = plot(k, c,
    label="Computed shooting path",
    xlabel="Capital kₜ",
    ylabel="Consumption cₜ",
    title="Saddle Path and Δk = 0 Curve",
    linewidth=2)
plot!(p1, kgrid, c_delta_k_zero,
    label="Δk = 0",
    linestyle=:dash,
    linewidth=2)
scatter!(p1, [kstar], [cstar],
    label="Steady state",
    markersize=4)
savefig(p1, "figures/saddle_path.png")

# =========================================================
# 8. PART (c): CAPITAL AND CONSUMPTION TIME SERIES
# =========================================================
p2 = plot(time, k,
    label="Capital kₜ",
    xlabel="Time t (periods)",
    ylabel="Value",
    title="Capital and Consumption Time Series",
    linewidth=2)
plot!(p2, time, c, label="Consumption cₜ", linewidth=2)
hline!(p2, [kstar], label="Steady-state capital", linestyle=:dash)
hline!(p2, [cstar], label="Steady-state consumption", linestyle=:dot)
savefig(p2, "figures/time_series.png")

# =========================================================
# 9. PART (d): CAPITAL POLICY PATH
# =========================================================
p3 = plot(k[1:T], k[2:T+1],
    label="Computed policy path",
    xlabel="Current capital kₜ",
    ylabel="Next-period capital kₜ₊₁",
    title="Capital Policy Path",
    linewidth=2,
    marker=:circle,
    markersize=2)
plot!(p3, kgrid, kgrid,
    label="45-degree line",
    linestyle=:dash)
savefig(p3, "figures/policy_function.png")

# =========================================================
# 10. PRINT RESULTS TO THE POWERSHELL WINDOW
# =========================================================
println("\n========== BEEM138 PROBLEM SET 4 ==========")
println("PART (a): Steady state")
@printf("k* = %.10f\n", kstar)
@printf("c* = %.10f\n", cstar)

println("\nPART (b): Shooting algorithm")
@printf("k0 = %.10f\n", k0)
@printf("c0 = %.12f\n", c0)
@printf("k_T = %.10f\n", k[end])
@printf("c_T = %.10f\n", c[end])

println("\nPART (c): Half-life")
@printf("Halfway capital target = %.10f\n", halfway_k)
if isnothing(half_life_period)
    println("Halfway point not reached within T periods.")
else
    println("Half-life = $half_life_period periods")
end

println("\nFiles saved:")
println("results/summary.txt")
println("results/time_series.csv")
println("results/policy_function.csv")
println("figures/saddle_path.png")
println("figures/time_series.png")
println("figures/policy_function.png")
println("============================================")

