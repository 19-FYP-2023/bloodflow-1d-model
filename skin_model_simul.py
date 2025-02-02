import sys
import numpy as np
import matplotlib
from configparser import SafeConfigParser
import matplotlib.pyplot as plt
from skinModel import model
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
    config.read(filename)

    order = config.getint('data', 'order')
    Nx = config.getint('data', 'Nx')
    Nt = config.getint('data', 'Nt')
    T0 = config.getfloat('data', 'T0')
    T = config.getfloat('data', 'T')
    L = [float(f) for f in config.get('data', 'L').split(',')]
    rc = config.getfloat('data', 'rc')
    qc = config.getfloat('data', 'qc')
    rho = config.getfloat('data', 'rho')
    mesh_locations = config.get('data', 'mesh_locations').split(',')
    names = config.get('data', 'names').split(',')
    locations = config.get('data', 'locations').split(',')

    return order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations,\
           names, locations



data_location = '/home/biyon/FYP/bloodflow-1d-model/patient_100_out/4cycles_last/data.cfg'

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
order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations, names, locations = read_output(data_location) 



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
m_d = 1613
k = (10 ** -6)
Qc = (328 * k) / (A*10)
print(Qc)
nd = 0.6875
na = 0.3125
tud = 0.39
tua = 1.78


j = 1 # (number of arteries = 3 => radial artery index = 1) | (number of arteries = 5 => radial artery index = 4) 
x = np.linspace(0, L[j], Nx+1)
print("length of atery: ", L[j])


D = [20,21.3,22.6]
for d in D:
    if d>L[j]:
        raise ValueError('Specified length is larger than atery length')
        
D = np.round(400*np.array(D)/L[j]).astype(int) # TODO: change 400 to Nx
print(D)

ppgs = []
for d in D:
    waveForm = dict()
    for i, name in enumerate(names):

            
        M = np.load('%s/%s_%i_M.npy' % (locations[i], name, j))
        print(M.shape)

        y = M[d, :]

        waveForm[name] = y

    #print(waveForm) #flow,area,pressure

    
    for key, val in waveForm.items():
        print(key, val.shape, type(val[0]))
        if (val.shape[0]!=Nt):
            raise ValueError('Waveforms Length is not consistant')


    """ 
    Simulating the flow rate effect for PPG
    """    
    diameter = 2*np.sqrt(waveForm["area"]/np.pi) # how the diameter change with time
 
    Q = (waveForm["flow"]/waveForm["area"])/100 # TODO: check if division by 100 is to convert cm to m


    y = np.sign(Q) * m_d * (np.sqrt(np.abs(Q / Qc)) / (1 + np.sqrt(np.abs(Q / Qc)))) # Eq,h is calculated here (paper: Quantification of the Phenomena Affecting Reflective Arterial Photoplethysmography)
    norm_y = y / np.max(y)

    norm_y_min = np.argmin(norm_y, axis=0)
    print("min p y : ", norm_y_min)

    norm_y_3 = np.tile(norm_y,3)

    delta = nd * np.exp(-time / tud) + na * np.exp(-time / tua) # related to equation (6) in the paper (paper: Quantification of the Phenomena Affecting Reflective Arterial Photoplethysmography)
    norm_delta = delta / np.sum(delta)

    ab_org = np.convolve(norm_y_3, norm_delta, 'same')
    ab_crop = ab_org[200:800] # TODO: experiment with the limit values

    ab_min = np.argmin(ab_crop, axis=0)
    print("min p y : ", ab_min)

    # check if the signal has been shifted forward after convolution
    if (ab_min<norm_y_min):
        raise ValueError("Min point cannot find")
    
    start_ab = 200 + ab_min - norm_y_min 
    # start_ab = 0 

    ab = ab_org[start_ab:start_ab+400] # TODO: replace 400 with Nt  

    # n = max(len(norm_y), len(norm_delta))
    # norm_y_centered = np.pad(norm_y, (n - len(norm_y), 0), mode='constant')
    # norm_delta_centered = np.pad(norm_delta, (n - len(norm_delta), 0), mode='constant')

    # # Convolve the centered signals
    # ab_centered = np.convolve(norm_y_centered, norm_delta_centered, 'same')

    # ab_croped = ab_centered[Nt:2*Nt] 
    print("min p ab : ", np.argmin(ab, axis=0))
     

    print('croped :',ab.shape)
    # plt.plot(norm_delta)
    # plt.plot(ab)
    # plt.plot(norm_y)

    # # plt.plot(time3, ab_centered)
    # # plt.plot(time, ab_croped)
    # plt.xlabel('t')
    # plt.ylabel('Q')
    # plt.title('abbsobtion from flow')
    # plt.grid(True)
    # plt.show()

    # print(diameter.shape)
    # print(waveForm["pressure"].shape)
    # print(ab.shape)

    """ 
    Skin simulation:

    - ab: Eq (in paper)

    """
    parameters = np.concatenate((diameter.reshape(Nt,1),waveForm["pressure"].reshape(Nt,1),ab.reshape(Nt,1)), axis=1)

    print(parameters.shape)

    nPhotonsCollected_values = np.zeros_like(time)

    # Loop over each time value
    for i in tqdm(range(len(time))):
        # Your existing code

        # print(parameters[i])
        nPhotonsCollected_values[i] = model(parameters[i])
        # print(nPhotonsCollected_values[i])

    # plt.plot(time, nPhotonsCollected_values, '-o')
    # plt.xlabel('t')
    # plt.ylabel('nPhotonsCollected')
    # plt.title('Simulation Results')
    # plt.grid(True)
    # plt.show()

    ppgs.append(nPhotonsCollected_values)
    plt.plot(time, nPhotonsCollected_values, label=f"D={d}")

plt.xlabel('Time')  # You may need to replace 'Time' with the appropriate label
plt.ylabel('Reflected light')  # You may need to replace 'Y' with the appropriate label
plt.title('Simulation Results')
plt.legend()
plt.show()