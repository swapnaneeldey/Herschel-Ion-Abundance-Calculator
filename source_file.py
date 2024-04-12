import csv

import numpy as np
from astropy.io import fits
"""This is the driver file which could read all the headers in the csv file and save its header contents 
as a callable function"""

def source_name():
    "returns all the source name as a list. egs [\"W3_88", "W3_122""...]"
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        source_name = []
        for row in reader:
            source_name.append(row[0])
    return source_name


def tau(source_name):
    """return the tau_9.7 for the called source. Saves a dict with the key is the name of the source with transition line and
    the value is the tau_9.7 of the source."""
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        tau_dict = {}
        for row in reader:
            tau_dict[row[0]] = float(row[1])
    return tau_dict[source_name]


def source_temperature(source_name):
    "returns the temperature for the called source. Saves a dict like the above function."
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        temperature_dict = {}
        for row in reader:
            temperature_dict[row[0]] = float(row[3])
    return temperature_dict[source_name]

def temperature_error(source_name):
    "returns the temperature error for the called source. Saves a dict like the above funtion."
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        temperature_error_dict = {}
        for row in reader:
            temperature_error_dict[row[0]] = float(row[4])
    return temperature_error_dict[source_name]
def ratio_eps(source_name):
    """returns the ratio of the N_e to N_p and also its error. Saves a dict like above but for value it
    saves the 1 + y and error in it."""
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        ratio_ep_dict = {}
        ratio_error_dict = {}
        for row in reader:
            ratio_ep_dict[row[0]] = 1 + float(row[6])
            ratio_error_dict[row[0]] = float(row[7])
    return ratio_ep_dict[source_name], ratio_error_dict[source_name]


def channel(source_name):
    """returns the good channel width of each transition line of each source which was inputed in the csv file manually.
    It saves a dict as well."""
    with open("/users/sdey/DATA/herschel.csv - Sheet2recent.csv", "r") as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(reader, None)
        channel_dict = {}
        for row in reader:
            first, second = row[8].split(":")
            channel_dict[row[0]] = (int(first),int(second))
    return channel_dict[source_name]




# ratio_ep("W3_88")

def bad_spaxels(source_name):
    """return the bad spaxel list, and it takes care of the bad spaxels due to 52 line. Using rms in the fits file and
    the channel changes for each source. Then it will give the bad spaxels of each source and its transition line in a
    form of a LIST. this list is then parsed through in the plotting codes."""

    position_file_name = source_name
    # get data
    # shape is (68,5,5) --> (channel, spaxel index, spaxel index)
    data = fits.getdata(f"/users/sdey/DATA/te_herschel/{position_file_name}_data.fits")
    # Exclude channels that cover the line profile (set to NaN)
    # This could be different for each source/transition; need to check.
    a,b = channel(source_name)
    data[a:b, :, :] = np.nan
    # exculde the band edges
    # data[65:,:,:] = np.nan
    # data[0:4,:,:] = np.nan
    # calculate the rms (5 x 5 array)
    rms = np.nanstd(data, axis=0)
    rms[[0, 4]] = rms[[4, 0]]
    rms[[1, 3]] = rms[[3, 1]]
    rms = rms.flatten()
    fit_params_file_path = f"/users/sdey/DATA/te_herschel/{position_file_name}_FitParams.txt"

    spaxel = [('4', '0'), ('4', '1'), ('4', '2'), ('4', '3'), ('4', '4'), ('3', '0'), ('3', '1'), ('3', '2'),
              ('3', '3'), ('3', '4'), ('2', '0'), ('2', '1'), ('2', '2'), ('2', '3'), ('2', '4'), ('1', '0'),
              ('1', '1'), ('1', '2'), ('1', '3'), ('1', '4'), ('0', '0'), ('0', '1'), ('0', '2'), ('0', '3'),
              ('0', '4')]

    with open(fit_params_file_path, "r") as fit_params_file:
        # Read each line in the FitParams file
        peak_intensity = []
        peak_intensity_dict = {}
        for line in fit_params_file:
            if "M1" in line:
                line_parts = line.split()
                line_position = tuple(line_parts[0:2])
                # Check if the line position matches any chosen spaxel
                if line_position in spaxel:
                    peak_intensity_dict[line_position] = float(line_parts[3])
    for spaxel_position in spaxel:
        peak_intensity.append(peak_intensity_dict[spaxel_position])

    #print(peak_intensity / rms)
    spaxel_to_skip = bad_spaxels_52(source_name)
    for n, (noise, signal) in enumerate(zip(rms, peak_intensity)):
        if signal / noise < 5 and spaxel[n] not in spaxel_to_skip:
            spaxel_to_skip.append(spaxel[n])
    return spaxel_to_skip

def bad_spaxels_52(source_name):
    """return the bad spaxel list for 52 transition line. Using rms in the fits file and the channel changes for each source. Then it will give
    the bad spaxels of each source and its transition line in a form of a LIST. this list is then parsed through in the plotting
    codes."""

    position_file_name = source_name.split("_")[0]
    # get data
    # shape is (68,5,5) --> (channel, spaxel index, spaxel index)
    data = fits.getdata(f"/users/sdey/DATA/te_herschel/{position_file_name}_52_data.fits")
    # Exclude channels that cover the line profile (set to NaN)
    # This could be different for each source/transition; need to check.
    a,b = channel(f"{position_file_name}_52")
    data[a:b, :, :] = np.nan
    # exculde the band edges
    # data[65:,:,:] = np.nan
    # data[0:4,:,:] = np.nan
    # calculate the rms (5 x 5 array)
    rms = np.nanstd(data, axis=0)
    rms[[0, 4]] = rms[[4, 0]]
    rms[[1, 3]] = rms[[3, 1]]
    rms = rms.flatten()
    fit_params_file_path = f"/users/sdey/DATA/te_herschel/{position_file_name}_52_FitParams.txt"

    spaxel = [('4', '0'), ('4', '1'), ('4', '2'), ('4', '3'), ('4', '4'), ('3', '0'), ('3', '1'), ('3', '2'),
              ('3', '3'), ('3', '4'), ('2', '0'), ('2', '1'), ('2', '2'), ('2', '3'), ('2', '4'), ('1', '0'),
              ('1', '1'), ('1', '2'), ('1', '3'), ('1', '4'), ('0', '0'), ('0', '1'), ('0', '2'), ('0', '3'),
              ('0', '4')]

    with open(fit_params_file_path, "r") as fit_params_file:
        # Read each line in the FitParams file
        peak_intensity = []
        peak_intensity_dict = {}
        for line in fit_params_file:
            if "M1" in line:
                line_parts = line.split()
                line_position = tuple(line_parts[0:2])
                # Check if the line position matches any chosen spaxel
                if line_position in spaxel:
                    peak_intensity_dict[line_position] = float(line_parts[3])
    for spaxel_position in spaxel:
        peak_intensity.append(peak_intensity_dict[spaxel_position])

    #print(peak_intensity / rms)
    spaxel_to_skip = []
    for n, (noise, signal) in enumerate(zip(rms, peak_intensity)):
        if signal / noise < 5:
            spaxel_to_skip.append(spaxel[n])
    return spaxel_to_skip
print(bad_spaxels("W3_122"))