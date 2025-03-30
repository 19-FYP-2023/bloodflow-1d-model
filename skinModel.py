from constants import *
import numpy as np
from scipy.stats import multivariate_normal
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from skin_model_const_params import SkinModelConstParams

def calc_mua(wavelength, S, B, W, F, M):
    global muadeoxy, muafat, muamel, muaoxy, muawater, musp, nmLIB
    if 'muadeoxy' not in globals():
        spectralLIB = np.load('spectralLIB.mat')
        muadeoxy = spectralLIB['muadeoxy']
        muafat = spectralLIB['muafat']
        muamel = spectralLIB['muamel']
        muaoxy = spectralLIB['muaoxy']
        muawater = spectralLIB['muawater']
        musp = spectralLIB['musp']
        nmLIB = spectralLIB['nmLIB']

    mua_deoxy = interp1d(nmLIB, muadeoxy)(wavelength)
    mua_fat = interp1d(nmLIB, muafat)(wavelength)
    mua_mel = interp1d(nmLIB, muamel)(wavelength)
    mua_oxy = interp1d(nmLIB, muaoxy)(wavelength)
    mua_water = interp1d(nmLIB, muawater)(wavelength)

    # Jacques "Optical properties of biological tissues: a review" eq. 12:
    mua = B * S * mua_oxy + B * (1 - S) * mua_deoxy + W * mua_water + F * mua_fat + M * mua_mel  # Equation 12 without the bilirubin and beta-Carotene terms
    return mua

def calc_mus(wavelength, aPrime, fRay, bMie, g):
    fMie = 1 - fRay
    musPrime = aPrime * (fRay * (wavelength / 500) ** (-4) + fMie * (wavelength / 500) ** (-bMie))  # Jacques "Optical properties of biological tissues: a review" eq. 2
    mus = musPrime / (1 - g)
    return mus


def mediaPropertiesFunc(parameters):
    class MediumProperties:
        def __init__(self, name, mua, mus, g, VHC, TC, E=None, A=None):
            self.name = name
            self.mua = mua
            self.mus = mus
            self.g = g
            self.VHC = VHC
            self.TC = TC
            self.E = E
            self.A = A

    def func_mua2(wavelength):
        B = 0  # Blood content
        S = 0.75  # Blood oxygen saturation
        W = 0.75  # Water content
        M = 0.03  # Melanin content
        F = 0  # Fat content
        if wavelength == 660:
            return 0.3442
        elif wavelength == 890:
            return 0.3184
        else:
            return calc_mua(wavelength, S, B, W, F, M)

    def func_mus2(wavelength):
        aPrime = 40  # musPrime at 500 nm
        fRay = 0  # Fraction of scattering due to Rayleigh scattering
        bMie = 1  # Scattering power for Mie scattering
        g = 0.9  # Scattering anisotropy
        if wavelength == 660:
            return 3.22
        elif wavelength == 890:
            return 224.7
        else:
            return calc_mus(wavelength, aPrime, fRay, bMie, g)

    def func_mua3(wavelength):
        B = 0.002  # Blood content
        S = 0.67  # Blood oxygen saturation
        W = 0.65  # Water content
        M = 0  # Melanin content
        F = 0  # Fat content
        if wavelength == 660:
            return 0.3162
        elif wavelength == 890:
            return 2459
        else:
            return calc_mua(wavelength, S, B, W, F, M)

    def func_mus3(wavelength):
        aPrime = 42.4  # musPrime at 500 nm
        fRay = 0.62  # Fraction of scattering due to Rayleigh scattering
        bMie = 1  # Scattering power for Mie scattering
        g = 0.9  # Scattering anisotropy
        if wavelength == 660:
            return 0.49
        elif wavelength == 890:
            return 116.7
        else:
            return calc_mus(wavelength, aPrime, fRay, bMie, g)

    def func_mua4(wavelength):
        B = 1  # Blood content
        S = 0.75  # Blood oxygen saturation
        W = 0.95  # Water content
        M = 0  # Melanin content
        F = 0  # Fat content
        if wavelength == 660:
            return 2.026
        elif wavelength == 890:
            return 6.32
        else:
            return calc_mua(wavelength, S, B, W, F, M)

    def func_mus4(wavelength):
        aPrime = 10  # musPrime at 500 nm
        fRay = 0  # Fraction of scattering due to Rayleigh scattering
        bMie = 1  # Scattering power for Mie scattering
        g = 0.9  # Scattering anisotropy
        if wavelength == 660:
            return 75.76
        elif wavelength == 890:
            return 56.18
        else:
            return calc_mus(wavelength, aPrime, fRay, bMie, g)

    def func_mua5(wavelength):
        B = 1  # Blood content
        S = 0.75  # Blood oxygen saturation
        W = 0.95  # Water content
        M = 0  # Melanin content
        F = 0  # Fat content
        if wavelength == 660 or wavelength == 890:
            return 0.8
        else:
            return calc_mua(wavelength, S, B, W, F, M)

    def func_mus5(wavelength):
        aPrime = 10  # musPrime at 500 nm
        fRay = 0  # Fraction of scattering due to Rayleigh scattering
        bMie = 1  # Scattering power for Mie scattering
        g = 0.9  # Scattering anisotropy
        if wavelength == 660 or wavelength == 890:
            return 230
        else:
            return calc_mus(wavelength, aPrime, fRay, bMie, g)

    def func_mua6(wavelength):
        B = 1  # Blood content
        S = 0.75  # Blood oxygen saturation
        W = 0.95  # Water content
        M = 0  # Melanin content
        F = 0  # Fat content
        if wavelength == 660:
            return 0.0001
        elif wavelength == 890:
            return 0.0217
        else:
            return calc_mua(wavelength, S, B, W, F, M)

    def func_mus6(wavelength):
        aPrime = 10  # musPrime at 500 nm
        fRay = 0  # Fraction of scattering due to Rayleigh scattering
        bMie = 1  # Scattering power for Mie scattering
        g = 0.9  # Scattering anisotropy
        if wavelength == 660:
            return 249.7
        elif wavelength == 890:
            return 189.8
        else:
            return calc_mus(wavelength, aPrime, fRay, bMie, g)

    mediaProperties = [
        MediumProperties('water', 0.00036, 10, 1.0, 4.19, 5.8e-3),
        MediumProperties('epidermis', func_mua2, func_mus2, 0.9, 3391*1.109e-3, 0.37e-2),
        MediumProperties('dermis', func_mua3, func_mus3, 0.9, 3391*1.109e-3, 0.37e-2),
        MediumProperties('blood', func_mua4, func_mus4, 0.9, 3617*1.050e-3, 0.52e-2, 422.5e3, 7.6e66),
        MediumProperties('vessel wall', func_mua5, func_mus5, 0.9, 3391*1.109e-3, 0.37e-2),
        MediumProperties('fat', func_mua6, func_mus6, 0.9, 3391*1.109e-3, 0.37e-2)
    ]

    return mediaProperties

# Helper functions calc_mua and calc_mus should be defined separately
# These functions are used inside the mediaPropertiesFunc but not provided in the given MATLAB code




def geometryDefinition(X, Y, Z, parameters):
    # Geometry parameters
    zsurf = 0.01
    epd_thick = 0.025
    dem_thick = 0.1
    vesselradius =  parameters[0] / 2
    vesseldepth = 0.375
    ves_wallthick = 0.02

    IDx = 0.4
    pL1 = 0.1
    pwd1 = 0.2

    IDxd = 0.3
    psep = 0.2
    pwd2 = 0.09
    frq = 660

    # Initialize matrix M with ones (background with water)
    M = np.ones_like(X) # FIX: replace X with Z 

    # Set regions based on Z coordinate
    M[Z > zsurf] = 2  # epidermis
    M[Z > zsurf + epd_thick] = 3  # dermis
    M[Z > zsurf + epd_thick + dem_thick] = 6  # fat

    # Define regions for vessel and vessel wall
    vessel_wall_indices = Y**2 + (Z - (zsurf + vesseldepth))**2 < (vesselradius + ves_wallthick)**2
    vessel_indices =  Y**2 + (Z - (zsurf + vesseldepth))**2 < vesselradius**2

    M[vessel_wall_indices] = 5  # vessel wall
    M[vessel_indices] = 4  # blood

    # Compute absorption coefficient A based on the defined regions
    mediaProperties = mediaPropertiesFunc([])
    vel_coeff = 10

    # print(mediaProperties[0].mua)

    A = np.ones_like(X) * mediaProperties[0].mua  # background with water (gel)
    A[Z > zsurf] = mediaProperties[1].mua(frq)  # epidermis
    A[Z > zsurf + epd_thick] = mediaProperties[2].mua(frq)  # dermis
    A[Z > zsurf + epd_thick + dem_thick] = mediaProperties[5].mua(frq)  # fat
    A[vessel_wall_indices] = mediaProperties[4].mua(frq)  # vessel wall
    A[vessel_indices] = mediaProperties[3].mua(frq) #+ vel_coeff * parameters[2]  # blood  # FIX: make the absorption coefficient reduce with Eq

    return M, A

# Assuming mediaPropertiesFunc is defined similarly to the MATLAB code
# You would need to implement mediaPropertiesFunc function for this code to work correctly
# You can reuse the previously defined mediaPropertiesFunc function if it's available

# Example usage:
# X, Y, Z = np.meshgrid(x_values, y_values, z_values)  # Define your grid coordinates
# parameters = [value1, value2, value3]  # Define your parameters
# M, A = geometryDefinition(X, Y, Z, parameters)

def model(parameters, skin_model_const_params, in_flux = 1):
    """
    MODEL Summary of this function goes here
    Detailed explanation goes here
    parameters
    1 - diameter
    """
    X = skin_model_const_params.get_meshX()
    Y = skin_model_const_params.get_meshY()
    Z = skin_model_const_params.get_meshZ()
    x = skin_model_const_params.get_x()
    y = skin_model_const_params.get_y()
    z = skin_model_const_params.get_z()

    # mx, my, mz = np.meshgrid(np.arange(1, 101))
    M, A = geometryDefinition(X, Y, Z, parameters)

    mean_light_pathway_z_points = skin_model_const_params.get_mean_light_path_z_vals()

    # setting the input and output flux values for the 1st x layer    
    xlayer_in_flux = in_flux
    xlayer_out_flux = None

    for i in range(len(x)):
        if skin_model_const_params.is_valid_mean_light_path_index(i):
            mu = [0, mean_light_pathway_z_points[i]] # [mean for y dir, mean for z dir]
            Sigma = [[RO * abs(mean_light_pathway_z_points[i]), 0], [0, RO * abs(mean_light_pathway_z_points[i])]] # [std for y dir, std for z dir]
            YM, ZM = np.meshgrid(y, z)
            YZ = np.column_stack((YM.ravel(), ZM.ravel()))
            norm = multivariate_normal.pdf(YZ, mu, Sigma).reshape(len(z), len(y))
            xlayer_voxels_mu = A[:, :, i]

            # calculating the flux going into each voxel in this x layer 
            xlayer_voxels_in_flux = xlayer_in_flux * norm * skin_model_const_params.get_voxel_yz_area()

            # calculating the flux coming out of each voxel in this x layer
            xlayer_voxels_out_flux = xlayer_voxels_in_flux * (1 - xlayer_voxels_mu * skin_model_const_params.get_ln10() * skin_model_const_params.get_mean_light_path_distances()[i])

            # calculating the total flux coming out of this x layer
            xlayer_out_flux = np.sum(xlayer_voxels_out_flux)

            # setting the input flux of the next layer as the output flux of this layer
            xlayer_in_flux = xlayer_out_flux

    # REVIEW: i think we have to show that the absSum we get here is equal to Edc + Ev + Eq in the paper
    return xlayer_out_flux