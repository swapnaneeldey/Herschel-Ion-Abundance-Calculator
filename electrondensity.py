import multiprocessing
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pyneb as pn
from operator import truediv

import scipy.stats
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
import effecient
import source_file

"""The only funtion getting called in this file is new_density and bad_density_spaxel function. 
 Other may be useful but depends how you use."""
o3 = pn.Atom("O", 3)

def density(source_name):
    """It gives out the electron density of every spaxel.
    Paramenters: Takes in the name of the source e.g. W3"""
    if source_name == "W51":
        return [4950]*25
    else:
        intensity_ratio = ratio(source_name)
        temp = source_file.source_temperature(f"{source_name}_88")
        electron_density_list = []
        for i in intensity_ratio:
            density = o3.getTemDen(i, tem=temp, wave1=520000, wave2=880000)
            electron_density_list.append(density)
        return electron_density_list


def ratio(source_name):
    corrected_flux_88 = effecient.corrected_flux(f"{source_name}_88")
    corrected_flux_52 = effecient.corrected_flux(f"{source_name}_52")
    intensity_ratio = list(map(truediv, corrected_flux_52, corrected_flux_88))
    return intensity_ratio

def bad_density_spaxels(source_name):
    bad_spaxels_88 = source_file.bad_spaxels(f"{source_name}_88")
    if source_name == "W51":
        bad_spaxels_52 = []
    else:
        bad_spaxels_52 = source_file.bad_spaxels(f"{source_name}_52")
    bad_spaxels_total = []
    for i in bad_spaxels_88:
        bad_spaxels_total.append(i)
    for j in bad_spaxels_52:
        if j not in bad_spaxels_total:
            bad_spaxels_total.append(j)
    return bad_spaxels_total

def formula(i,j):
    density = o3.getTemDen(i, tem=j, wave1=520000, wave2=880000)
    if not np.isnan(density) and np.isfinite(density):
        return density
    else:
        return 1

def new_density(source_name):
    start = time()
    print("starting electron Density calculation")
    corrected_flux_88 = np.array(effecient.corrected_flux(f"{source_name}_88"))
    corrected_flux_52 = np.array(effecient.corrected_flux(f"{source_name}_52"))
    corrected_flux_88_error = np.array(effecient.corrected_flux_error(f"{source_name}_88"))
    corrected_flux_52_error = np.array(effecient.corrected_flux_error(f"{source_name}_52"))
    intensity_ratio = corrected_flux_52 / corrected_flux_88
    intensity_ratio_error = intensity_ratio * (np.sqrt(
        (corrected_flux_88_error / corrected_flux_88) ** 2 + (corrected_flux_52_error / corrected_flux_52) ** 2))

    # Define the Gaussian function for curve fitting
    def gaussian(x, mu, sigma, amplitude):
        return amplitude * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    mean_electron_density = []
    monte_carlo_density_error = []
    for i in range(0, len(intensity_ratio_error)):
        s = time()
        flux = np.abs(np.random.normal(intensity_ratio[i], intensity_ratio_error[i], 4000))
        temp = source_file.source_temperature(f"{source_name}_88")
        temp_error = source_file.temperature_error(f"{source_name}_88")
        temperature = np.random.normal(temp, temp_error, 4000)
        with multiprocessing.Pool(16) as pool:
            electron_density_list = list(pool.starmap(formula, zip(flux,temperature)))
        electron_density_list = [density for density in electron_density_list if density != 1]
        data = np.array(electron_density_list)
        kde = KernelDensity(kernel="gaussian", bandwidth=75).fit(data[:, np.newaxis])
        x_grid = np.linspace(np.min(electron_density_list), np.max(electron_density_list), 1100)[:,np.newaxis]
        log_dens = kde.score_samples(x_grid)
        peak, _ = find_peaks(np.exp(log_dens))
        height = np.exp(log_dens[peak])
        print(x_grid[np.argmax(np.exp(log_dens))])


        # Perform curve fitting on the KDE data
        #p0 = [float(x_grid[np.argmax(np.exp(log_dens))]), float(scipy.stats.median_abs_deviation(data)), np.max(x_grid)]
        # Initial guess for curve fitting parameters
        p0 = [np.median(data), scipy.stats.median_abs_deviation(data), np.max(np.exp(log_dens))]
        print(p0)
        params, cov_matrix = curve_fit(gaussian, x_grid.flatten(), np.exp(log_dens), p0=p0)
        fitted_curve = gaussian(x_grid, *params)

        mean_electron_density.append(params[0])
        monte_carlo_density_error.append(params[1])
        e = time()
        print(e - s)
        end = time()
        print(end - start)
        print(np.array(monte_carlo_density_error)/np.array(mean_electron_density))
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax1.hist(electron_density_list, density=True, bins=85)
    ax1.set_xlabel("electron_density")
    ax1.set_ylabel("count")
    ax1.plot(x_grid, np.exp(log_dens), label="KDE gaussian", color="red")
    ax1.plot(x_grid, fitted_curve, label='Fitted Gaussian', color='green')
    ax1.legend()
    print(np.mean(data))
    ax1.set_title(np.mean(np.array(electron_density_list)))
    ax2.hist(temperature, bins=85)
    ax2.set_xlabel("temperature")
    ax3.hist(flux, bins=85)
    ax3.set_xlabel("flux ratio")
    #plt.show()
    return np.array(mean_electron_density), np.array(monte_carlo_density_error)

#print(new_density("W3"))

"""peak, _ = find_peaks(np.exp(log_dens))
print(int(peak))
height = np.exp(log_dens[peak])/2
left_index = np.argmin(np.abs((np.exp(log_dens[:int(peak)]) - height)))
right_index = np.argmin(np.abs(np.exp(log_dens[int(peak):]) - height))
fwhm = x_grid[right_index + int(peak)]-x_grid[left_index]
sigma = fwhm/2.35
center = x_grid[peak]
print(center)
print(fwhm)
print(sigma)"""
"""fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
ax1.hist(electron_density_list, density = True, bins = 55)
ax1.set_xlabel("electron_density")
ax1.set_ylabel("count")
ax1.plot(x_grid, np.exp(log_dens), label = "KDE gaussian", color = "red")
ax1.plot(x_grid, fitted_curve, label='Fitted Gaussian', color='green')
ax1.legend()
print(np.mean(data))
ax1.set_title(np.mean(np.array(electron_density_list)))
ax2.hist(temperature, bins=55)
ax2.set_xlabel("temperature")
ax3.hist(flux, bins=55)
ax3.set_xlabel("flux ratio")
plt.show()"""