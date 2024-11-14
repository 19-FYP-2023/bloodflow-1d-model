# import os
print("Test Image")
# print("Current Location :", os.getcwd())
# print(os.listdir())

import dolfinx

packages_to_check = ["FiniteElement", "FunctionSpace", "FunctionSpace_float64", "Expression", "Function", "TestFunctions", "TestFunction", "DirichletBC", "solve"]
available_package_list = dir(dolfinx.fem)

for package in packages_to_check:
    if package in available_package_list:
        print(f"{package} is available")
    else:
        print(f"[ERROR] {package} is not available")
