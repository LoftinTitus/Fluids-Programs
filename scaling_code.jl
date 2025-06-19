# Scaling Code
# Code for calculating basic scaling laws on pulsatile flows with assmetric step of ~50% of the diameter

# Constant integers
# Two different radii
r_real = 0.01/2
r_scaled = 0.025/2

# Frequency of the flow (70 BPM)
f = 1.167

# Dynamic viscoity
mu = 3.5e-3

# Density of the fluid
rho = 1060

# Calculated for Womersley number
alpha_real = r_real * sqrt((2 * pi * f * mu) / rho)
alpha_scaled = r_scaled * sqrt((2 * pi * f * mu) / rho)

Q = [2.45, 1.63, 2.37, 0.7, 1.4, 1.54, 1.25, 1.89, 1.49] # List of experimental flow rates (L/min)

# This function takes an array of experimental flow rates and turns them into their respective avg velcoites
function Q_to_u_scaled(Q, r_scaled)
    u_scaled = Float64[]
    A = pi * r_scaled^2

    for q in Q
        # Convert from L/min to m³/s: L/min * (1m³/1000L) * (1min/60s)
        q_m3_s = (q * (1/1000)) * (1/60)
        # Calculate velocity: Q = A * u, so u = Q/A
        push!(u_scaled, q_m3_s/A)
    end

    return u_scaled
end 

# Calculate u_scaled using the function
u_scaled = Q_to_u_scaled(Q, r_scaled)

# Initialize arrays for storing results
tau_array = Float64[]
tau_reduced_array = Float64[]
v_equiv_array = Float64[]
tau_reduced_real_array = Float64[]
tau_p_real_array = Float64[]
tau_min_real_array = Float64[]

# Calculated for Reynolds number
function Re_scaled(u_scaled, r_scaled, mu, rho)
    Re_values = Float64[]
    for velocity in u_scaled
        push!(Re_values, (rho * velocity * 2 * r_scaled) / mu)
    end
    return Re_values
end

# This function calculates the scaled shear stresses
function min_tau(Re_scaled, u_scaled)
    # constants for the scaling law
    a = 10.476
    b = 2.127e-2
    c = -9.776e-6

    min_tau_values = Float64[]
    for (re, u) in zip(Re_scaled, u_scaled)
        # Calculate the shear stress using the scaling law

        tau_poly = a * re * b + c * re^2
        tau_min = tau_poly * (2 * mu * u) / r_scaled
        push!(min_tau_values, tau_min)
    end
    return min_tau_values
end

# This function caluclates the posieelle shear stress
function posieulle_tau(r_scaled, mu, u_scaled)
    posieulle_tau_values = Float64[]
    for velocity in u_scaled
        tau_posieulle = (4* mu * velocity) / r_scaled
        push!(posieulle_tau_values, tau_posieulle)
    end
    return posieulle_tau_values
end

# Calculate Reynolds numbers
Re_values = Re_scaled(u_scaled, r_scaled, mu, rho)

# Calculate minimum tau values
min_tau_values = min_tau(Re_values, u_scaled)

# Calculate Poiseuille tau values
posieulle_tau_values = posieulle_tau(r_scaled, mu, u_scaled)

# Debug output
println("u_scaled: ", u_scaled)
println("Re_values: ", Re_values)
println("min_tau_values: ", min_tau_values)
println("posieulle_tau_values: ", posieulle_tau_values)

# This function calculates the reduced shear stress and shear stress in pa using NL solve
using NLsolve

solutions = []

# Loop through each set of values
for (i, u) in enumerate(u_scaled)
    t_min = min_tau_values[i]
    t_p = posieulle_tau_values[i]
    
    # Constant Declarations
    global r = r_scaled
    global mu = mu

    # System of equations
    function system!(F, x)
        # List of my unknowns
        tau = x[1]
        tau_reduced = x[2]

        # My two equations
        # Add small epsilon to avoid division by zero
        eps = 1e-12
        denominator = t_min - t_p
        if abs(denominator) < eps
            denominator = eps
        end
        
        F[1] = tau_reduced - (tau - t_p)/denominator
        F[2] = tau - ((2 * mu * u) / r) * (tau_reduced * (t_min - 2) + 2)
    end

    # Initial guess for [tau, tau_reduced]
    x0 = [1.0, 0.5]
    sol = nlsolve(system!, x0)

    tau_sol = sol.zero[1]
    tau_reduced_sol = sol.zero[2]

    # Push into separate arrays
    push!(tau_array, tau_sol)
    push!(tau_reduced_array, tau_reduced_sol)

end

# This function is going to solve for the equivalent velocity after scaling down
for (i, v) in enumerate(tau_array)
    tau = tau_array[i]

    global r_real = r_real
    global mu = mu
    global rho = rho

    # constants for the scaling law
    a = 10.476
    b = 2.127e-2
    c = -9.776e-6

    # Reynolds number coefficent ie: Re  = d * u
    d = (rho * r_real * 2) / mu

    function system!(F, x)

        # List of my unknowns
        v_equiv = x[1]
        tau_reduced_real = x[2]
        tau_p_real = x[3]
        tau_min_real = x[4]

        # Add small epsilon to avoid division by zero
        eps = 1e-12
        denominator = tau_min_real - tau_p_real
        if abs(denominator) < eps
            denominator = eps
        end

        # My four equations
        F[1] = tau - ((2 * mu * v_equiv) / r_real) * (tau_reduced_real * (tau_min_real - 2) + 2)
        F[2] = tau_reduced_real - (tau - tau_p_real)/denominator
        F[3] = tau_p_real - (4 * mu * v_equiv) / r_real
        F[4] = tau_min_real - (a + d * v_equiv * b + c * d^2 * v_equiv^2)
    end

    # Initial guess for [v_equiv, tau_reduced_real, tau_p_real, tau_min_real]
    x0 = [1.0, 0.5, 0.5, 0.5]
    sol = nlsolve(system!, x0)
    v_equiv_sol = sol.zero[1]
    tau_reduced_real_sol = sol.zero[2]
    tau_p_real_sol = sol.zero[3]
    tau_min_real_sol = sol.zero[4]

    # Push into separate arrays
    push!(v_equiv_array, v_equiv_sol)
    push!(tau_reduced_real_array, tau_reduced_real_sol)
    push!(tau_p_real_array, tau_p_real_sol)
    push!(tau_min_real_array, tau_min_real_sol)
end

# Print our array of velocities scaled down to realistic size
println("Scaled velocities: ", v_equiv_array)

# Convert velocitys into flow rates in L/min
function u_to_Q_real(v_equiv_array, r_real)
    Q_real = Float64[]
    A = pi * r_real^2

    for v in v_equiv_array

        push!(Q_real, A * v * 1000 * 60) # Convert m³/s to L/min
    end

    return Q_real
end

Q_scaled = u_to_Q_real(v_equiv_array, r_real)

println("Scaled flow rates:", Q_scaled)