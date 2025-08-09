import scipy.io
import numpy as np
import pandas as pd


mat_data = scipy.io.loadmat('/home/dumindu/modeling/pwdb_data.mat')


#get patient number from system arguments
import sys
if len(sys.argv) < 2:
    print("Usage: python createCardiacOutput.py <patient_number>")
    sys.exit(1)

patient = int(sys.argv[1]) - 1  # Convert to zero-based index
if patient < 0 or patient >= mat_data['data'][0,0]['haemods'].shape[1]:
    print(f"Invalid patient number: {patient + 1}. Please provide a valid patient number.")
    sys.exit(1)

age  = mat_data['data'][0,0]['haemods'][0,patient]['age'][0,0]
print(f"Processing patient {patient}...in age {age}")

flowRateDate = mat_data['data'][0,0]['waves'][0,0]['U_AorticRoot'][0,patient].flatten()

# old_min, old_max = flowRateDate.min(), flowRateDate.max()
# new_min, new_max = 2, 25

# # # Scale the data
# scaledFlowRateDate = ((flowRateDate) * (new_max - new_min)) / (old_max)
#scaledFlowRateDate = flowRateDate

# duplicate array and concatenate to itself n times
n = 1
scaledFlowRateDate = np.tile(flowRateDate, n)

fs=mat_data['data'][0,0]['waves'][0,0]['fs'][0][0]

# Generate time column
num_samples = scaledFlowRateDate.shape[0]
time = np.arange(num_samples) / fs

# Create a DataFrame
data = pd.DataFrame({
    "Time (s)": time,
    "Value": scaledFlowRateDate
})


# Save to CSV with a more appropriate filename
output_filename = f"data/Original_aortic/input.csv"
data.to_csv(output_filename, index=False, header=False)

print(f"Data has been saved for patient :'{patient}'")
