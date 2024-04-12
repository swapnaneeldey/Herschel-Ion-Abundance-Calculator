import numpy as np

import source_file as sf

"This file gets the flux of the sources. The main part of the code is making a big dictionary which has the" \
"correct_integral, integral, np.exp(tau_88), true_intensity, correct_integral_error, true_intensity_error" \
"for all the sources. " \
"Then a smaller functions call this big dictionary to get the specific item from a specific source."

def calculate_total_intensity(peak_intensity, sigma_lambda, error_sigma_lambda, error_peak_intensity, transition_line,
                              position_file_name, tau_9_7):
    """
    Calculates the total intensity based on peak intensity and sigma lambda values.
    """
    integral = []
    if transition_line == 52:
        # correction factor for weak 52 line detection by herschel.
        correction_factor = 273.62 / 15.13
        for intensity, sigma in zip(peak_intensity, sigma_lambda):
            # Calculate the intensity value using the given formula
            value = ((3 * 10 ** 8) / (transition_line * 10 ** (-6)) ** 2) * 1.06 * intensity * 2.35 * sigma * 10 ** (
                -6) * 10 ** -23 * correction_factor
            integral.append(value)

        integral_error = integral * np.sqrt(
            (np.divide(error_peak_intensity, peak_intensity)) ** 2 + (np.divide(error_sigma_lambda, sigma_lambda)) ** 2
            + (3.73 / 273.62) ** 2 + (4.081067 / 1513) ** 2)

        return (extinction(transition_line, integral, position_file_name, tau_9_7, integral_error))
    else:
        for intensity, sigma in zip(peak_intensity, sigma_lambda):
            # Calculate the intensity value using the given formula
            value = ((3 * 10 ** 8) / (transition_line * 10 ** (-6)) ** 2) * 1.06 * intensity * 2.35 * sigma * 10 ** (
                -6) * 10 ** -23
            integral.append(value)

        integral_error = np.sqrt(
            (np.divide(error_peak_intensity, peak_intensity)) ** 2 + (np.divide(error_sigma_lambda, sigma_lambda)) ** 2)

        return (extinction(transition_line, integral, position_file_name, tau_9_7, integral_error))


def extinction(transition_line, integral, position_file_name, tau_9_7, integral_error):
    """Function that calculates the corrected flux of every spaxel by taking into account the extinctions of it in
    every transition. Here observed intensity is the sum of the integral"""
    observed_intensity = sum(integral)
    # from literature
    # formula to get the optical depth of the dust at a specific wavelength
    tau_52 = 0.054 * tau_9_7
    tau_57 = 0.044 * tau_9_7
    tau_88 = 0.019 * tau_9_7
    tau_122 = 0.0098 * tau_9_7
    dict_flux = {}
    # makes a dictionary for every transtion for every source as its key and for the values it has one list of
    # correct fluxes, another list for observed fluxes, an e^tau_wavelength value, sum of the corrected fluxes.
    if transition_line == 88:
        correct_integral = [np.exp(tau_88) * i for i in integral]
        correct_integral_error = [np.sqrt(i**2 + (tau_88 * 0.20)**2 + 0.15**2)*j for i,j in zip(integral_error, correct_integral)]
        true_intensity = observed_intensity * np.exp(tau_88)
        true_intensity_error = np.sqrt(sum(np.array(correct_integral_error) ** 2))
        #dict is being made
        dict_flux[position_file_name] = correct_integral, integral, np.exp(tau_88), true_intensity, \
            correct_integral_error, true_intensity_error
        return dict_flux
    elif transition_line == 52:
        correct_integral = [np.exp(tau_52) * i for i in integral]
        correct_integral_error = [np.sqrt(i**2 + (tau_52 * 0.20)**2 + 0.15**2) * j for i,j in zip(integral_error, correct_integral)]
        true_intensity = observed_intensity * np.exp(tau_52)
        true_intensity_error = np.sqrt(sum(np.array(correct_integral_error) ** 2))
        dict_flux[position_file_name] = correct_integral, integral, np.exp(tau_52), true_intensity, \
            correct_integral_error, true_intensity_error
        return dict_flux
    elif transition_line == 57:
        correct_integral = [np.exp(tau_57) * i for i in integral]
        correct_integral_error = [np.sqrt(i**2 + (tau_57 * 0.20)**2 + 0.15**2) * j for i,j in zip(integral_error, correct_integral)]
        true_intensity = observed_intensity * np.exp(tau_57)
        true_intensity_error = np.sqrt(sum(np.array(correct_integral_error) ** 2))
        dict_flux[position_file_name] = correct_integral, integral, np.exp(tau_57), true_intensity, \
            correct_integral_error, true_intensity_error
        return dict_flux
    else:
        true_intensity = observed_intensity * np.exp(tau_122)
        correct_integral = [np.exp(tau_122) * i for i in integral]
        correct_integral_error = [np.sqrt(i**2 + (tau_122 * 0.20)**2 + 0.15**2) * j for i,j in zip(integral_error, correct_integral)]
        true_intensity_error = np.sqrt(sum(np.array(correct_integral_error) ** 2))
        dict_flux[position_file_name] = correct_integral, integral, np.exp(tau_122), true_intensity, \
            correct_integral_error, true_intensity_error
        return dict_flux



def correct_flux(position_file_name):
    """Function to parse through fitsparam txt file to get the peak intensity
    width of the gaussian fit."""

    file_name = "Spaxel"
    # with open(file_name, "r") as position_file:
    # Read each line in the position file
    # for position in position_file:
    # for position_file_name in sf.source_name():
    # position_values = position.split(" , ")
    # position_file_name = position_values[0]
    tau_9_7 = sf.tau(position_file_name)
    # print(tau_9_7)
    # Extract the transition line number from the file name
    transition_line = int(position_file_name.split("_")[1])

    # Initialize lists to store peak intensity and sigma lambda values and errors
    peak_intensity = []
    sigma_lambda = []
    error_peak_intensity = []
    error_sigma_lambda = []

    spaxel = [('4', '0'), ('4', '1'), ('4', '2'), ('4', '3'), ('4', '4'), ('3', '0'), ('3', '1'), ('3', '2'),
              ('3', '3'), ('3', '4'), ('2', '0'), ('2', '1'), ('2', '2'), ('2', '3'), ('2', '4'), ('1', '0'),
              ('1', '1'), ('1', '2'), ('1', '3'), ('1', '4'), ('0', '0'), ('0', '1'), ('0', '2'), ('0', '3'),
              ('0', '4')]

    # Iterate over the spaxel pairs in the spaxel file
    """for i in range(1, len(position_values)):
        x, y = position_values[i][0], position_values[i][2]
        #spaxel.append((x, y))"""
    peak_intensity_dict = {}
    sigma_lambda_dict = {}
    error_peak_intensity_dict = {}
    error_sigma_lambda_dict = {}
    # Open the corresponding FitParams file
    fit_params_file_path = f"/users/sdey/DATA/te_herschel/{position_file_name}_FitParams.txt"
    with open(fit_params_file_path, "r") as fit_params_file:
        # Read each line in the FitParams file
        for line in fit_params_file:
            if "M1" in line:
                line_parts = line.split()
                line_position = tuple(line_parts[0:2])
                # Check if the line position matches any chosen spaxel
                if line_position in spaxel:
                    # Extract the peak intensity and sigma lambda values
                    try:
                        peak_intensity_dict[line_position] = float(line_parts[3])
                        error_peak_intensity_dict[line_position] = float(line_parts[4])
                        sigma_lambda_dict[line_position] = (float(line_parts[7]))
                        error_sigma_lambda_dict[line_position] = (float(line_parts[8]))
                    except ValueError:
                        # print("skipping")
                        spaxel.remove(line_position)
                    if not spaxel:
                        break  # Exit the loop if all unique positions have been found
    # print(spaxel)
    for spaxel_position in spaxel:
        peak_intensity.append(peak_intensity_dict[spaxel_position])
        sigma_lambda.append(sigma_lambda_dict[spaxel_position])
        error_peak_intensity.append(error_peak_intensity_dict[spaxel_position])
        error_sigma_lambda.append(error_sigma_lambda_dict[spaxel_position])
    # print("spaxels taken", spaxel)
    # Calculate and print the total intensity for the position
    return (calculate_total_intensity(peak_intensity, sigma_lambda, error_sigma_lambda, error_peak_intensity,
                                      transition_line, position_file_name, tau_9_7))


#important
list_dic = []
for name in ["W3_88", "W3_52", "W3_122", "W3_57" ]:
    # makes one list with all the dictionary of every source transition made above.
    list_dic.append(correct_flux(name))

print(list_dic)

#functions that get called in different codes.
def corrected_flux(source_name):
    """A function that can be called to get the corrected flux list of a particular transition of a source."""
    for u in list_dic:
        if source_name in u:
            return u[source_name][0]

def not_corrected_flux(source_name):
    """A function that can be called to get the observed flux list of a particular transition of a source."""
    for u in list_dic:
        if source_name in u:
            return u[source_name][1]


def corrected_flux_error(source_name):
    "Return a list of error per spaxel for the corrected flux."
    for u in list_dic:
        if source_name in u:
            return u[source_name][4]

