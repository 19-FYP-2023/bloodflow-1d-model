import scipy.io
import numpy as np
import pandas as pd

# Load MATLAB file
mat_data = scipy.io.loadmat('/Users/biyon/Documents/MySpace/Repo/pwdb_data.mat')


patient = 100

print(mat_data.keys())
print(mat_data['data'].shape)
print(mat_data['data'][0,0].dtype)
print(mat_data['data'][0,0]['waves'].shape)
print(mat_data['data'][0,0]['waves'][0,0].dtype)
print(mat_data['data'][0,0]['waves'][0,0]['U_AorticRoot'].shape)
print(mat_data['data'][0,0]['waves'][0,0]['U_AorticRoot'][0,patient].shape)

flowRateDate = mat_data['data'][0,0]['waves'][0,0]['U_Brachial'][0,patient].flatten()

old_min, old_max = flowRateDate.min(), flowRateDate.max()
new_min, new_max = 2, 25

# Scale the data
scaledFlowRateDate = ((flowRateDate - old_min) * (new_max - new_min)) / (old_max - old_min)


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
output_filename = f"data/patient_{patient}_flow_rate.csv"
data.to_csv(output_filename, index=False, header=False)

print(f"Data has been saved to '{output_filename}'")






# patient = 40
# cross_section = np.array(mat_data['data']['waves']['A_Radial'][0][0][patient][0])
# velocity = 50 * np.array(mat_data['data']['waves']['U_Radial'][0][0][patient][0]) + 1
# pressure = np.array(mat_data['data']['waves']['P_Radial'][0][0][patient][0])
# PPG = np.array(mat_data['data']['waves']['PPG_Radial'][0][0][patient][0])

# Diameter = np.sqrt((4/np.pi) * cross_section)

# # Assuming you have the sampling frequency (fs) available
# fs = mat_data['data']['waves']['fs'][0][0][0][0]

# # Calculate time vector based on the number of samples and sampling frequency
# num_samples = len(pressure)
# time = np.arange(num_samples) / fs

# # Plotting the waveforms
# fig, axs = plt.subplots(4, 1, figsize=(10, 12))

# axs[0].plot(time, Diameter, 'b', linewidth=2)
# axs[0].set_title('Cross Section Diameter Waveform')
# axs[0].set_xlabel('Time (s)')
# axs[0].set_ylabel('Cross Section Diameter (m)')
# axs[0].grid(True)

# axs[1].plot(time, velocity, 'g', linewidth=2)
# axs[1].set_title('Velocity Waveform')
# axs[1].set_xlabel('Time (s)')
# axs[1].set_ylabel('Velocity (cm/s)')
# axs[1].grid(True)

# axs[2].plot(time, pressure, 'r', linewidth=2)
# axs[2].set_title('Pressure Waveform')
# axs[2].set_xlabel('Time (s)')
# axs[2].set_ylabel('Pressure (mmHg)')
# axs[2].grid(True)

# new_num_samples = 145
# new_time = np.linspace(time[0], time[-1], new_num_samples)
# new_velocity = np.interp(new_time, time, velocity)

# # Plotting the resampled waveforms
# axs[1].plot(new_time, new_velocity, 'g', linewidth=2)
# axs[1].set_title('Velocity Waveform (Resampled)')
# axs[1].set_xlabel('Time (s)')
# axs[1].set_ylabel('Velocity (cm/s)')
# axs[1].grid(True)

# # Save data to CSV file
# with open('example_inlet.csv', 'w') as file:
#     for i in range(len(new_time)):
#         file.write(f'{new_time[i]}, {new_velocity[i]}\n')

# print('CSV file successfully generated: example_inlet.csv')
# plt.tight_layout()
# plt.show()
