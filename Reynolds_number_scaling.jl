# Reynolds number scaling
# This program is for scaling via keeping Reynolds number constants

# Include original scaling code to use Re function
include("scaling_code.jl")

# Print array of already calulated Re numbers
println("Scaled Reynolds number:", Re_values)

# Function to solve for the real fluid velocities after keeping Reynolds number constant
function solve_for_velocity(r_scaled, r_real, u_scaled)
    real_velocity_values = Float64[]  # Declare the array
    
    # Because mu, and rho are constant, we can use the following relation:
    # u_real = u_scaled * (r_scaled / r_real)
    for u in u_scaled
        u_real = u * (r_scaled / r_real)  # Use individual element u, not entire array
        push!(real_velocity_values, u_real)
    end

    return real_velocity_values  # Fixed typo
end

real_velocity_values = solve_for_velocity(r_scaled, r_real, u_scaled)

println("Real fluid velocities after scaling:", real_velocity_values)

# Convert velocites to flowrates
real_Q_values = u_to_Q_real(real_velocity_values, r_real)

println("Real flow rates after scaling:", real_Q_values)