import concurrent
import multiprocessing
from time import time

import numpy as np
from operator import truediv

from astropy.io import fits

import effecient
import freefreecalc
import source_file
import emissivity as em
import source_file as sf

"""This code helps in getting the ionic abundance and the error of the source for each of it transitions.  This code has 3 functions which
will produce results for all the results together and store it in the form of a dictionary, which is later called by smaller function to get the 
abundance/abundance_error of a specific source."""


def error(value, flux_error, corrected_flux, ratio_error, ratio_ep, temp_error, temperature, emissivity,
          emissivity_error, freefree_error, free_free):
    log_abundance_error = value * np.sqrt(
        (flux_error / corrected_flux) ** 2 + (ratio_error / ratio_ep) ** 2 + (
                temp_error * 0.35 / temperature) ** 2 + (emissivity_error / emissivity) ** 2 + (
                    freefree_error / free_free) ** 2)

    return log_abundance_error


def formula(corrected_flux, constant, temperature, frequency, ratio_ep, free_free, emissivity, flux_error, ratio_error,
            temp_error, emissivity_error, freefree_error):
    value = (corrected_flux * constant * (temperature ** -0.35) * (frequency ** -0.1) * ratio_ep) / (
            free_free * emissivity)
    log_abundance_error = error(value, flux_error, corrected_flux, ratio_error, ratio_ep, temp_error, temperature,
                                emissivity,
                                emissivity_error, freefree_error, free_free)
    # return float("{:.2e}".format(np.log10(value) + 12)), log_abundance_error
    return float("{:.2e}".format(value)), float("{:.2e}".format(log_abundance_error))


abundance_dict = {}


def abundance():
    """makes the big dictionary with all the information, i.e, {source1:(abundance, abundance error), source2:..}"""
    s = time()
    for name in ["W3_88", "W3_122", "W3_52",
                 "W3_57"]:  # replace the list with sf.source_name() for automating it for all possible sources.
        print(f"Starting {name}")
        # print(dict)
        # try:
        pixel_val = freefreecalc.freefree(name)
        # except FileNotFoundError:
        # continue
        corrected_flux = effecient.corrected_flux(name)
        # print("f", corrected_flux[12])
        flux_error = effecient.corrected_flux_error(name)
        # print(flux_error[12])
        ratio_ep, ratio_error = sf.ratio_eps(name)
        # print(ratio_ep, ratio_error)
        frequency = em.frequency(name)
        temperature = sf.source_temperature(name) / 10000
        # print(temperature)
        temp_error = sf.temperature_error(name) / 10000
        # print(temp_error)
        emissivity, emissivity_error = em.emissivity_range(name)
        # print(emissivity, emissivity_error)
        constant = 3.485 * (10 ** (-16))
        freefree_error = 15 * 10 ** -3
        # print("fe", flux_error[12] / corrected_flux[12], ratio_error / ratio_ep, temp_error * 0.35 / temperature,
        # emissivity_error[12] / emissivity[12])
        abundance_list = []
        abundance_error_list = []
        for i, (key, value) in enumerate(pixel_val.items()):
            free_free = value
            abundances, abundance_error = formula(corrected_flux[i], constant, temperature, frequency, ratio_ep,
                                                  free_free,
                                                  emissivity[i], flux_error[i], ratio_error, temp_error,
                                                  emissivity_error[i], freefree_error)
            abundance_list.append(abundances)
            abundance_error_list.append(abundance_error)
        abundance_dict[name] = abundance_list, abundance_error_list
    e = time()
    print("time taken: ", e - s)

    return abundance_list, abundance_error_list


# fills up the dict #important
abundance()


"""def multiabundance():
    if __name__ == "__main__":
        names = ["W3_52", "W3_57"]

        with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
            futures = {executor.submit(abundance, name): name for name in names}


            for future in concurrent.futures.as_completed(futures):
                name = futures[future]
                result = future.result()
                abundance_dict[name] = result

        return abundance_dict

multiabundance()"""


def source_abundance(source_name):
    return [x for x in abundance_dict[source_name][0]]


def source_abundance_error(source_name):
    return abundance_dict[source_name][1]


# print(source_abundance_error("W3_88")[12])

def abundance_mean(source_name):
    abundance = source_abundance(source_name)
    abundance_error = source_abundance_error(source_name)
    total_weight = 0
    weighted_sum = 0

    for value, error in zip(abundance, abundance_error):
        weight = 1 / (error ** 2)
        total_weight += weight
        weighted_sum += value * weight

    error = np.sqrt(1 / total_weight)

    return weighted_sum / total_weight, error


# print(abundance_mean("W3_88"))


#fllowing is not useful right now
def nitrogen_ratio(source):
    return (list(
        map(truediv, abundance_dict[f"{source}_57"][0], abundance_dict[f"{source}_122"][0]))), skip_spaxels_nitrogen(
        source)


def nitroxy_ratio(source):
    return list(
        map(truediv, abundance_dict[f"{source}_57"][0], abundance_dict[f"{source}_88"][0])), skip_spaxels_nitroxy(
        source)


def skip_spaxels_nitroxy(source_name):
    bad_spaxels_88 = source_file.bad_spaxels(f"{source_name}_88")
    bad_spaxels_57 = source_file.bad_spaxels(f"{source_name}_57")
    bad_spaxels_total = []
    for i in bad_spaxels_88:
        bad_spaxels_total.append(i)
    for j in bad_spaxels_57:
        if j not in bad_spaxels_total:
            bad_spaxels_total.append(j)
    return bad_spaxels_total


def skip_spaxels_nitrogen(source_name):
    bad_spaxels_122 = source_file.bad_spaxels(f"{source_name}_122")
    bad_spaxels_57 = source_file.bad_spaxels(f"{source_name}_57")
    bad_spaxels_total = []
    for i in bad_spaxels_122:
        bad_spaxels_total.append(i)
    for j in bad_spaxels_57:
        if j not in bad_spaxels_total:
            bad_spaxels_total.append(j)
    return bad_spaxels_total
