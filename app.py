# import os
print("Test Image")
# print("Current Location :", os.getcwd())
# print(os.listdir())

import dolfinx

packages_to_check = ["FiniteElement", "FunctionSpace", "Expression", "Function", "TestFunctions", "DirichletBC", "solve"]
available_package_list = dir(dolfinx.fem)

for package in packages_to_check:
    if package in available_package_list:
        print(f"{package} is available")
    else:
        print(f"{package} is not available")
