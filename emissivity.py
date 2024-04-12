import multiprocessing
from glob import glob
from time import time

import pyneb as pn
import numpy as np
import scipy
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from sklearn.neighbors import KernelDensity

import electrondensity
import source_file

"""This code caluclates the emissivity of the source for each spaxel.Only 2 main functions are being used, frequency and 
 emissivity range."""



#print(multiprocessing.cpu_count())


def frequency(name):
    """get the radio frequency from vla radio fits file and divides it by 5GHz"""
    source = name.split("_")[0]
    hdu = fits.open(f"/users/sdey/DATA/{source}_radio.fits")
    header = hdu[0].header
    return header["CRVAL3"] / (5 * 10 ** 9)


def emissivity(name):
    o3 = pn.Atom("O", 3)
    n2 = pn.Atom("N", 2)
    n3 = pn.Atom("N", 3)
    source = name.split("_")[0]
    transition_line = int(name.split("_")[1])
    density = electrondensity.density(source)
    temp = source_file.source_temperature(f"{source}_88")
    if transition_line == 122:
        em = n2.getEmissivity(tem=temp, den=density, wave=1220000)
    elif transition_line == 57:
        em = n3.getEmissivity(tem=temp, den=density, wave=570000)
    else:
        em = o3.getEmissivity(tem=temp, den=density, wave=transition_line * 10000)
    return list(em)


o3 = pn.Atom("O", 3)
n2 = pn.Atom("N", 2)
n3 = pn.Atom("N", 3)

def transition_lines(name):
    global transition_line
    transition_line = (int(name.split("_")[1]))
    return transition_line


def formula(d, t):
    """For multiprocessing"""
    if transition_line == 122:
        em = n2.getEmissivity(tem=t, den=d, wave=1220000)
    elif transition_line == 57:
        em = n3.getEmissivity(tem=t, den=d, wave=573237.70)
    else:
        em = o3.getEmissivity(tem=t, den=d, wave=transition_line * 10000)
    return em

def emissivity_range(name):
    """Calculates the emissivity and emissivity error list for each spaxel using monte carlo method. It also prints out a plot for
    visualizing the KDE and Gaussian fit of the KDE."""
    start = time()
    source = name.split("_")[0]
    #making the transition line a global variable so it can be used in the multiprocessing.
    transition_lines(name)

    #getting the density and density error from the electrondensity.py.
    density, density_error = electrondensity.new_density(source) #if you go to the midspaxel.txt file i have copy pasted the electro density values of W3 there to save time.
    temp = source_file.source_temperature(f"{source}_88")
    temp_error = source_file.temperature_error(f"{source}_88")
    mean_emissivity = []
    monte_carlo_emissivity_error = []
    print("starting emissivity calculation")
    for i in range(0, len(density)):
        s = time()
        den_normal = np.abs(np.random.normal(density[i], density_error[i], 3000))
        temp_normal = np.random.normal(temp, temp_error, 3000)
        with multiprocessing.Pool(16) as pool:
            em_list = list(pool.starmap(formula, zip(den_normal, temp_normal)))
        #dividing by 10**-22 so its more efficient in calculating
        em_list = np.array(em_list) / 10 ** -22
        data = np.array(em_list)
        print(data)
        print(len(em_list))
        x_grid = np.linspace(np.min(em_list), np.max(em_list), 1100)[:, np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=0.05).fit(data[:, np.newaxis])
        log_em = kde.score_samples(x_grid)
        print(x_grid)

        def gaussian(x, mu, sigma, amplitude):
            return amplitude * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

        # Perform curve fitting on the KDE data
        p0 = [np.mean(data), np.std(data), np.max(data)]  # Initial guess for curve fitting parameters
        print(p0)
        params, _ = curve_fit(gaussian, x_grid.flatten(), np.exp(log_em), p0=p0)
        print(params)
        # Generate the fitted curve using the optimized parameters
        fitted_curve = gaussian(x_grid, *params)

        sigma = params[1]
        mean = params[0]
        mean_emissivity.append(mean * 10 ** -22)
        monte_carlo_emissivity_error.append(sigma * 10 ** -22)
        e = time()
        plt.show()
        print(e - s)
        print(np.array(monte_carlo_emissivity_error) / np.array(mean_emissivity))
    #Will Plot the fitting of the last spaxel in iteration.
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax1.hist(em_list, density=True, bins=75)
    ax1.set_xlabel("emissivity")
    ax1.set_ylabel("count")
    ax1.plot(x_grid, np.exp(log_em), label="KDE gaussian", color="red")
    ax1.plot(x_grid, fitted_curve, label='Fitted Gaussian', color='green')
    ax1.legend()
    print(np.mean(data))
    ax1.set_title(np.mean(np.array(em_list))*10**-22)
    ax2.hist(temp_normal, bins=85)
    ax2.set_xlabel("temperature")
    ax3.hist(den_normal, bins=85)
    ax3.set_xlabel("flux ratio")
    #plt.show()
    end = time()
    print(end - start)
    return mean_emissivity, monte_carlo_emissivity_error




print(emissivity_range("W3_88"))


