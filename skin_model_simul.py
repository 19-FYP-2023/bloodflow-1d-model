from configparser import SafeConfigParser
from constants import *
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import signal
from skinModel import model
from skin_model_const_params import SkinModelConstParams
import sys
from tqdm import tqdm

def read_output(filename):
    """
    Read data file generated in the output folder.

    Arguments
    ---------
    rc : string
        Data file name

    Returns
    -------
    return : tuple
        Tuple of all parameters stored in the file
    """
    config = SafeConfigParser()
    config.read(filename+'4cycles_last/data.cfg')

    no_of_arteries = config.getint('data', 'no_of_arteries')
    Nx = config.getint('data', 'nx')
    Nt = config.getint('data', 'nt')
    T0 = config.getfloat('data', 't0')
    T = config.getfloat('data', 't')
    L = [float(f) for f in config.get('data', 'l').split(',')]
    rc = config.getfloat('data', 'rc')
    qc = config.getfloat('data', 'qc')
    rho = config.getfloat('data', 'rho')
    mesh_locations = config.get('data', 'mesh_locations').split(',')
    names = config.get('data', 'names').split(',')
    locations = config.get('data', 'locations').split(',')

    for i,loc in enumerate(mesh_locations):
        mesh_locations[i] = filename+loc 

    for i,loc in enumerate(locations):
        locations[i] = filename+loc 


    return no_of_arteries, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations,\
           names, locations


if len(sys.argv) != 2:
    print("[ERROR] Incorrect invocation of the script")
    print("Usage: python skin_model_simul.py <patient data folder name>")
    sys.exit(1)

# making a folder with the name of the patient if not available
patient_folder_name = sys.argv[1]
base = f"outputs/{patient_folder_name}/"
result_folder = f"results/{patient_folder_name}/"

if not os.path.exists(base):
    print(f"[ERROR] Patient data does not exist in {base}")
    sys.exit(1)

if not os.path.exists(result_folder):
    print(f"[LOG] creating {result_folder}")
    os.mkdir(f"{result_folder}")
else:
    print(f"[LOG] cleaning existing {result_folder}")
    os.rmdir(f"{result_folder}")
    os.mkdir(f"{result_folder}")

""" 
used variables: Nx, Nt, T0, T, L, names, locations 

- Nx: number of divisions from the artery length to get the delta_x value
- Nt: number of divisions from the (T0 -> T) time range to get the delta_t value
- T0: start time of the simulation
- T: end time of the simulation
- L: artery lengths array (if there are 3 arteries then there will be a list of three values)
- names: 'area', 'flow', 'pressure'
- locations: locations of where the blood flow simulations outputs are for area, flow, and pressure
"""
no_of_arteries, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations, names, locations = read_output(base) 

# REVIEW: both Nx and L is read from data. So the length referred to by a unit of Nx can vary. Is that an issue?  
print(f"T0 = {T0}, T = {T}, Nt = {Nt}")
sampling_freq = Nt/(T-T0)
nyquist_freq = sampling_freq/2
norm_cutoff_freq = CUTOFF_FREQ/nyquist_freq

print(f"Sampling Frequency = {sampling_freq} Hz")
print(f"Nyquist Frequency = {nyquist_freq} Hz")

time = np.linspace(T0, T, Nt)
time3 = np.linspace(T0, T0+3*(T-T0), 3*Nt)
fs = Nt/(T-T0)

""" 
skin model parameters: 

A:
m_d:
k:
Qc:
nd:
na:
tud:
tua:
"""
A = np.pi * (( 1.75 / 1000) ** 2) 
m_d = 1613 # REVIEW: check whether this value is correct
k = (10 ** -6)
Qc = (328 * k) / (A*10) # REVIEW: check whether this value is correct
print(Qc)
nd = 0.6875
na = 0.3125
tud = 0.39
tua = 1.78


j = 4 # (number of arteries = 3 => radial artery index = 1) | (number of arteries = 5 => radial artery index = 4) 
x = np.linspace(0, L[j], Nx+1)
print("length of atery: ", L[j])


D = [22]
for d in D:
    if d>L[j]:
        raise ValueError('Specified length is larger than atery length')
        
D = np.round(Nx*np.array(D)/L[j]).astype(int)
print(D)

ppgs = []
for d in D:
    waveForm = dict()
    for i, name in enumerate(names):

            
        M = np.load('%s/%s_%i_M.npy' % (locations[i], name, j))
        print(M.shape)

        y = M[d, :]

        waveForm[name] = y
    
    for key, val in waveForm.items():
        print(key, val.shape, type(val[0]))
        if (val.shape[0]!=Nt):
            raise ValueError('Waveforms Length is not consistant')
    
    #plot area,flowrate and pressure of wave form in three subplots
    fig, axs = plt.subplots(3)
    fig.suptitle('Waveform')
    axs[0].plot(time, waveForm["area"])
    axs[0].set_title('Area')
    axs[1].plot(time, waveForm["flow"])
    axs[1].set_title('Flow')
    axs[2].plot(time, waveForm["pressure"])
    axs[2].set_title('Pressure')
    fig.savefig(f"{result_folder}/input_waveforms_d({d}).png")

    """ 
    Simulating the flow rate effect for PPG
    """    
    diameter = 2*np.sqrt(waveForm["area"]/np.pi) # how the diameter change with time

    """ 
    Skin simulation:

    - ab: Eq (in paper)

    """
    # parameters has the following time varying quantities at a specific location:
    #   - diameter of the vessel
    #   - pressure inside the vessel
    #   - effect of blood flow for rPPG
    parameters = np.concatenate((diameter.reshape(Nt,1),waveForm["pressure"].reshape(Nt,1)), axis=1)

    print(parameters.shape)

    nPhotonsCollected_values = np.zeros_like(time)
    eql_vals = np.zeros_like(time)
    eqh_vals = np.zeros_like(time)

    skin_model_const_params = SkinModelConstParams(GEO_XLIMITS, GEO_YLIMITS, GEO_ZLIMITS, NUM_X_TICKS, NUM_Y_TICKS, NUM_Z_TICKS, TX_RX_DISTANCE, MEAN_PENETRATION)

    # Loop over each time value
    for i in tqdm(range(len(time))):
        # propagate light through the skin model and get the remaining amount arriving at the receiver
        nPhotonsCollected_values[i] = model(parameters[i], skin_model_const_params)

    np.save(f"{result_folder}/ppg_d({d}).npy", nPhotonsCollected_values)

    # perform low-pass filtering with a butterworth filter to remove high frequency noise
    b, a = signal.butter(FILTER_ORDER, norm_cutoff_freq, btype='low')
    filtered_inv_signal = signal.filtfilt(b, a, nPhotonsCollected_values)

    filtered_ppg_signal = 1 - filtered_inv_signal
    ppgs.append(filtered_ppg_signal)

    np.save(f"{result_folder}/filtered_ppg_signal({d}).npy", filtered_ppg_signal)

    fig_2, ax_2 = plt.subplots()
    fig_2.suptitle(f"PPG_d({d})")
    ax_2.plot(time, filtered_ppg_signal)
    fig_2.savefig(f"{result_folder}/ppg_d({d}).png")

plt.xlabel('Time')  # You may need to replace 'Time' with the appropriate label
plt.ylabel('Reflected light')  # You may need to replace 'Y' with the appropriate label
plt.title('Simulation Results')
plt.legend()
plt.savefig(f"{result_folder}/final.png")