import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.wcs import WCS
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord, FK4, FK5, Galactic
from numpy import cos
import source_file

"""Automated footprint of herschel overlay on the radio image for each source in epoch J2000."""

for name in source_file.source_name():
    file = name
    source = file.split("_")[0]
    transition_line = int(file.split("_")[1])
    if transition_line == 88:
        try:
            hdu = fits.open(f"/users/sdey/DATA/{source}_radio.fits")
        except FileNotFoundError:
            continue
        header = hdu[0].header
        data = hdu[0].data
        data = data[0, 0]
        wcs = WCS(header)
        # drop freq and pol axes
        wcs = wcs.dropaxis(3)
        wcs = wcs.dropaxis(2)
        #
        # plot image
        #
        # zoom in
        ra_b1950, dec_b1950 = wcs.wcs.crval

        # Convert the B1950 coordinates to J2000
        coord_b1950 = SkyCoord(ra=ra_b1950 * u.deg, dec=dec_b1950 * u.deg, frame='fk4', equinox='B1950.0')
        coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))

        # Update the WCS to J2000
        wcs.wcs.crval = coord_j2000.ra.deg, coord_j2000.dec.deg
        #
        fig, ax = plt.subplots()
        ax = plt.subplot(projection=wcs)
        # ax = plt.subplot(projection=wcscut, slices=('x','y',0,0))
        im = ax.imshow(data)
        # cn = ax.contour(data, levels=[0.2], colors='k')
        cbar = fig.colorbar(im)
        cbar.set_label(r'Flux Density [Jy beam$^{-1}$]')
        #
        # Display KAO beam only for W3 Source.
        # Rudolph et al. (2006, ApJS, 162, 346)
        # W3A: B1950: RA = 02 21 56.3; Dec = 61 52 47
        ra = (2.0 + 21.0 / 60.0 + 56.3 / 3600.0) * 15
        dec = 61.0 + 52.0 / 60.0 + 47.0 / 3600.0
        # print(ra, dec)
        beam = 40.0 / 3600.0
        coord_b1950 = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='fk4', equinox='B1950.0')
        coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))
        ra = coord_j2000.ra.deg
        dec = coord_j2000.dec.deg
        # print(ra, dec)
        r = SphericalCircle((ra, dec * u.deg), beam / 2.0 * u.degree, edgecolor='red',
                            transform=ax.get_transform('fk4'), facecolor='none')

        ax.add_patch(r)
        #
        # Display Herschel footprint
        spaxel = 9.4 / 3600.0

        # transforming the coordinates of the W3 into 1950
        data_cube = f"/users/sdey/DATA/te_herschel/{file}_dec.fits"

        cube = fits.open(data_cube)
        cube_header = cube[0].header
        cube_data = cube[1].data
        dec = cube_data[0]
        # print(dec)
        dec[[0, 4]] = dec[[4, 0]]
        dec[[1, 3]] = dec[[3, 1]]
        decData = dec.flatten()
        dec_j2000 = np.round(decData, 4)
        # print("dec_j2000", dec_j2000)

        data_cube = f"/users/sdey/DATA/te_herschel/{file}_ra.fits"

        cube = fits.open(data_cube)
        cube_header = cube[0].header
        cube_data = cube[1].data
        ra = cube_data[0]
        # print(ra)
        ra[[0, 4]] = ra[[4, 0]]
        ra[[1, 3]] = ra[[3, 1]]
        raData = ra.flatten()
        ra_j2000 = np.round(raData, 4)
        # print("ra_j2000", ra_j2000)

        increment = 0
        x_val = []
        y_val = []
        # loop through to overlay the spaxels onto the radio image.
        for i, (ra_val, dec_val) in enumerate(zip(raData, decData)):
            x, y = wcs.world_to_pixel_values(raData[i], decData[i])
            x_val.append(x)
            y_val.append(y)
            r = SphericalCircle((raData[i] * u.deg, decData[i] * u.deg), spaxel / 2.0 * u.degree, edgecolor='white',
                                transform=ax.get_transform('fk4'), facecolor='none')
            # print(raData[i], decData[i])
            ax.add_patch(r)
            ax.text(x, y, str(4 - increment) + ", " + str(i % 5), color="white", ha='center', va='center')
            # print(str(4 - increment) + ", " + str(i % 5))
            if i % 5 == 4:
                increment += 1

        if transition_line == 88:
            plt.title(f"{source}")

        elif transition_line == 52:
            plt.title(f"{source}_O[III]_{str(transition_line)}")

        elif transition_line == 57:
            plt.title(f"{source}_N[III]_{str(transition_line)}")

        else:
            plt.title(f"{source}_N[III]_{str(transition_line)}")

        # save image
        # fig.tight_layout()
        plt.xlabel("RA (J2000)")
        plt.ylabel("Dec (J2000)")
        plt.grid(True)
        plt.gca().set_aspect('equal')
        # print(min(raData), max(raData))
        # print(min(decData), max(decData))
        ax.set_xlim(min(x_val) - 15, max(x_val) + 15)
        ax.set_ylim(min(y_val) - 15, max(y_val) + 15)
        plt.savefig(f"/users/sdey/DATA/j2000_plots/{source}")
        plt.show()
