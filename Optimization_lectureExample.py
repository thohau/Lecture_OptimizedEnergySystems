from scipy.optimize import linprog


"""A A company has to buy 100 MWh electrical energy at the market, they have two offers:
Material Costs in €/MWh CO2 emission in t/MWh
Coal 50 €/MWh and 0.9 t C02/MWh
Gas 90 €/MWh and 0.4 t C02/MWh
The CO2 emission for the 100 MWh must be limited to maximal 60 t. Formulate a linear
optimization problem to find the lowest price meeting all requirements."""
from scipy.optimize import linprog

# Goal: Minimization of costs
# Objective function coefficients (cost per MWh for Coal and Gas)
# We aim to minimize the total cost.
c = [50, 90]

# Inequality constraints (A_ub * x <= b_ub):
# We need to ensure that the CO2 emissions do not exceed 60 tons.
# This is represented by the constraint 0.9 * x1 + 0.4 * x2 <= 60
A_ub = [[0.9, 0.4]]
b_ub = [60]

# Equality constraints (A_eq * x = b_eq):
# The total amount of energy needed is 100 MWh.
# This is represented by the constraint x1 + x2 = 100
A_eq = [[1, 1]]
b_eq = [100]

# Bounds for variables:
# Both x1 (energy from coal) and x2 (energy from gas) should be non-negative.
# This is represented by the bounds (0, None) which means x1 >= 0 and x2 >= 0
x0_bounds = (0, None)
x1_bounds = (0, None)

# Calling the solver:
# We use the linprog function from scipy.optimize to solve the linear programming problem.
# We pass in the objective function coefficients, inequality constraints, equality constraints, and variable bounds.
result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[x0_bounds, x1_bounds], method='highs')

# Outputting the results:
# Check if the solver found a successful solution
if result.success:
    print("Optimal Solution:")
    # The optimal amount of energy to be obtained from coal
    print("x1 (MWh from Coal):", result.x[0])
    # The optimal amount of energy to be obtained from gas
    print("x2 (MWh from Gas):", result.x[1])
    # Since the problem is a minimization problem, we directly use the result.
    print("Minimum cost in Euro:", result.fun)
else:
    # If the solver did not find a solution, we print an appropriate message
    print("No solution was found.")