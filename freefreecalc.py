import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK4
import astropy.units as u
from scipy.ndimage import map_coordinates

def freefree(source_name):
    """Calculates the bremsstrahlung radiation at the center of each spaxel for a source.
    Returns a dictionary where the keys are the spaxel number and value is the free free radiation value
    e.g. {("4", "0"): 0.234}"""
    # read the VLA image.
    file = source_name
    source = file.split("_")[0]
    hdu = fits.open(f"/users/sdey/DATA/{source}_radio.fits")
    header = hdu[0].header
    data = hdu[0].data
    image_data = data[0, 0]
    wcs = WCS(header)[0, 0]

    # getting the coordinate of herschel pointing at the source.
    data_cube = f"/users/sdey/DATA/te_herschel/{file}_dec.fits"

    cube = fits.open(data_cube)
    cube_header = cube[0].header
    cube_data = cube[1].data
    dec = cube_data[0]
    # print(dec)
    dec[[0, 4]] = dec[[4, 0]]
    dec[[1, 3]] = dec[[3, 1]]
    dec = dec.flatten()
    dec_j2000 = np.round(dec, 4)

    data_cube = f"/users/sdey/DATA/te_herschel/{file}_ra.fits"

    cube = fits.open(data_cube)
    cube_header = cube[0].header
    cube_data = cube[1].data
    ra = cube_data[0]
    # print(ra)
    ra[[0, 4]] = ra[[4, 0]]
    ra[[1, 3]] = ra[[3, 1]]
    ra = ra.flatten()
    ra_j2000 = np.round(ra, 4)

    coord_j2000 = SkyCoord(ra=ra_j2000 * u.deg, dec=dec_j2000 * u.deg, frame="fk5", equinox="J2000")
    coord_b1950 = coord_j2000.transform_to(FK4(equinox="B1950"))
    ra = coord_b1950.ra.deg
    dec = coord_b1950.dec.deg
    #print(ra)
    #0print(dec)

    increment = 0
    dict = {}

    # calculating the free free

    for i, (ra_val, dec_val) in enumerate(zip(ra, dec)):
        # getting the center of the spaxel of herschel in terms of pixel value of the radio image.
        x, y = wcs.world_to_pixel_values(ra_val, dec_val)
        #print(x, y)
        #print(ra_val,dec_val)
        #print(str(4 - increment) + ", " + str(i % 5))
        key = (str(4 - increment) + ", " + str(i % 5))
        value = (float(x), float(y))
        dict[key] = value
        if i % 5 == 4:
            increment += 1
    #print(dict)
    pixel_val_dict = {}
    for key, value in dict.items():
        x, y = value
        # Perform bilinear interpolation to get the pixel value
        # Prepare the coordinate arrays
        coords = np.array([[y], [x]])

        # Perform bilinear interpolation to get the pixel value
        pixel_value = map_coordinates(image_data, coords, order=0, mode='nearest')[0]
        #print(pixel_value)
        pixel_val_dict[key] = pixel_value
    return pixel_val_dict

