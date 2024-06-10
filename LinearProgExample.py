# Sample from
# https://realpython.com/linear-programming-python/
#
# Maximize z = x + 2 y 
# Constraints (colors according to tutorial)
#  2x +  y <= 20  (red) 
# -4x + 5y <= 10  (blue)
#  -x + 2y >= -2  (yellow, multiply with -1)
#
# Install scipy in advance:
# python3 -m pip install scipy

# Load linprog from scipy
from scipy.optimize import linprog

# objective function (coefficients of linear function
obj = [-1, -2]
#      --  --
#       |   |- Coefficient for y
#       |----- Coefficient for x

# Constraint matrix
lhs_ineq = [[ 2,  1],  # Red constraint left side
            [-4,  5],  # Blue constraint left side
            [ 1, -2]]  # Yellow constraint left side

# Constratin rhs
rhs_ineq = [20,  # Red constraint right side
            10,  # Blue constraint right side
             2]  # Yellow constraint right side

# Equality constraint
lhs_eq = [[-1, 5]]  # Green constraint left side
rhs_eq = [15]       # Green constraint right side

# Limits for x and y
bnd = [(0, float("inf")),  # Bounds of x
       (0, float("inf"))]  # Bounds of y
	  
# Run Simples (linprog) solver	  
# method="revised simplex" is deprecated, use "highs" instead.
opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq,
              A_eq=lhs_eq, b_eq=rhs_eq, bounds=bnd,
              method="highs")
			  
# Show results			  
print(opt)
