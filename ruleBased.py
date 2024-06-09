# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np

# Data
hourlyPrice = [21.75, 30.29, 13.49, 30.14, 21.00, 17.32, 17.91, 18.30, 28.77, 14.97, 32.24, 21.56, 40.66, 41.70, 19.89, 16.66, 23.40, 22.98, 30.60, 30.09, 31.76, 13.21, 28.36, 30.53]
COPHeatPump = [2.0, 2.41, 2.82, 3.21, 3.60, 3.96, 4.31, 4.62, 4.90, 5.15, 5.35, 5.52, 5.64, 5.70, 5.72, 5.68, 5.59, 5.44, 5.24, 4.97, 4.65, 4.28, 3.84, 3.35]

heatDemandBuilding = 50  # kWh per day
powerHeatpumpElectrical = 3  # kW

# Compute efficiency factor (COP / price)
efficiency_factor = [COPHeatPump[i] / hourlyPrice[i] for i in range(24)]

# Sort hours by efficiency factor in descending order
sorted_hours = sorted(range(24), key=lambda i: efficiency_factor[i], reverse=True)

# Initialize operation times to zero
operation_times = [0] * 24
total_heat_produced = 0

# Allocate operation times based on efficiency
for hour in sorted_hours:
    if total_heat_produced >= heatDemandBuilding:
        break
    operation_times[hour] = 1  # Run heat pump at full capacity
    total_heat_produced += powerHeatpumpElectrical * COPHeatPump[hour]

# Calculate total cost
total_cost = sum(hourlyPrice[i] * operation_times[i] * powerHeatpumpElectrical for i in range(24))

# Print results
total_cost = round(total_cost / 100, 2)
print("\nTotal cost", total_cost)

# Data for visualization
hours = range(1, 25)

# Thermal energy generation of the heat pump as hourly array
thermal_energy = [operation_times[i] * powerHeatpumpElectrical * COPHeatPump[i] for i in range(24)]

# Output of thermally generated energy
print("Set point of the heat pump")
for i, val in enumerate(operation_times):
    print(f"Hour {i+1}: {val*100:.2f} %")

# Output of thermally generated energy
print("Thermal energy generation of the heat pump (kWh per hour):")
for i, val in enumerate(thermal_energy):
    print(f"Hour {i+1}: {val:.2f} kWh")

sum_thermal_energy = round(sum(thermal_energy), 2)
print("Sum of generated thermal energy:")
print(sum_thermal_energy)

# Electrical power consumption of the heat pump
plt.figure(figsize=(12, 6))
plt.bar(hours, thermal_energy, label='Thermal power of heat pump in kW', color='green', alpha=0.7)

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
