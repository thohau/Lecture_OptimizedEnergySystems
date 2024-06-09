# Import necessary libraries
from scipy.optimize import minimize  # For optimization
import matplotlib.pyplot as plt      # For plotting
import numpy as np                    # For numerical operations

# Data
hourlyPrice = [21.75, 30.29, 13.49, 30.14, 21.00, 17.32, 17.91, 18.30, 28.77, 14.97, 32.24, 21.56, 40.66, 41.70, 19.89, 16.66, 23.40, 22.98, 30.60, 30.09, 31.76, 13.21, 28.36, 30.53]
#hourlyPrice = [21.75, 3.29, 13.49, 30.14, 21.00, 17.32, 17.91, 18.30, 28.77, 14.97, 32.24, 21.56, 40.66, 41.70, 19.89, 16.66, 23.40, 22.98, 30.60, 30.09, 3.76, 13.21, 28.36, 30.53]

COPHeatPump = [2.0, 2.41, 2.82, 3.21, 3.60, 3.96, 4.31, 4.62, 4.90, 5.15, 5.35, 5.52, 5.64, 5.70, 5.72, 5.68, 5.59, 5.44, 5.24, 4.97, 4.65, 4.28, 3.84, 3.35]

heatDemandBuilding = 50  # kWh per day
powerHeatpumpElectrical = 3  # kW

# Objective and constraint functions
def objective_function(x):
    total_cost = 0
    for i in range(len(x)):
        total_cost += hourlyPrice[i] * x[i] * powerHeatpumpElectrical
    return total_cost

def constraint(x):
    total_heat_produced = sum(powerHeatpumpElectrical * COPHeatPump[i] * x[i] for i in range(len(x)))
    return heatDemandBuilding - total_heat_produced

# Initial guess for optimization
initial_guess = [0] * 24

# Bounds for variables (0 to 1 to limit operation time between 0 and 100%)
bounds = [(0, 1)] * 24

# Optimization
result = minimize(objective_function, initial_guess, bounds=bounds, constraints={'type': 'eq', 'fun': constraint})

# Total cost
totalCost = round(result.fun/100, 2)

print("\nTotal cost", totalCost)

# Data for visualization
hours = range(1, 25)

# Optimized operation time of the heat pump
optimized_hours = result.x

# Thermal energy generation of the heat pump as hourly array
thermal_energy = [optimized_hours[i] * powerHeatpumpElectrical * COPHeatPump[i] for i in range(len(hours))]


# Output of thermally generated energy
print("Set point of the heat pump")
for i, val in enumerate(optimized_hours):
    print(f"Hour {i+1}: {val*100:.2f} %")


# Output of thermally generated energy
print("Thermal energy generation of the heat pump (kWh per hour):")
for i, val in enumerate(thermal_energy):
    print(f"Hour {i+1}: {val:.2f} kWh")

sum_ThermalEnergy = round(sum(thermal_energy),2)

print("Sum of generated thermal energy:")
print(sum_ThermalEnergy)

# Electrical power consumption of the heat pump
plt.figure(figsize=(12, 6))



plt.bar(hours, thermal_energy, 
        label='Thermal power of heat pump in kW', color='green', alpha=0.7)

# Second y-axis on the right for electricity price
plt.twinx()

# Electricity price as line with secondary y-axis
plt.plot(hours, hourlyPrice, label='Electricity price in Cent/kWh', color='red', linestyle='--')

# Hourly COP
plt.plot(hours, COPHeatPump, label='COP heat Pump', color='grey', linestyle='--')

# First y-axis on the left for energy quantity
plt.ylabel('Electricity price')

# Plot properties
plt.title('Optimized operation time of the heat pump and electricity price')
plt.xlabel('Hour of the day')
plt.xticks(range(1, 25))
plt.grid(True)

# Legends
lines, labels = plt.gca().get_legend_handles_labels()
lines2, labels2 = plt.gcf().get_axes()[0].get_legend_handles_labels()
plt.legend(lines + lines2, labels + labels2, loc='upper right')

# Show plot
plt.show()
