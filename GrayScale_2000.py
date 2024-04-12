import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import SkyCoord, FK4, FK5
import matplotlib.cm as cm

import effecient
import electrondensity
import source_file as sf

"This code is automated in a different python file called Automated_Grayscale.py"
def grayscale(file, Type):
    """Makes the grayscale of each transisition of each source.
     Parameters: File: essentialy the source name with its transition line.
     Type: A string saying the type of grayscale <Abundance, Flux, Density>"""

    # reading the radio image.
    file = file
    grayscale = Type
    source = file.split("_")[0]
    transition_line = int(file.split("_")[1])
    hdu = fits.open(f"/users/sdey/DATA/{source}_radio.fits")
    header = hdu[0].header
    data = hdu[0].data
    data = data[0, 0]
    wcs = WCS(header)
    # drop freq and pol axes
    wcs = wcs.dropaxis(3)
    wcs = wcs.dropaxis(2)

    ra_b1950, dec_b1950 = wcs.wcs.crval
    # Convert the B1950 coordinates to J2000
    coord_b1950 = SkyCoord(ra=ra_b1950 * u.deg, dec=dec_b1950 * u.deg, frame='fk4', equinox='B1950.0')
    coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))

    # Update the WCS to J2000
    wcs.wcs.crval = coord_j2000.ra.deg, coord_j2000.dec.deg
    ax = plt.subplot(projection=wcs)
    # ax = plt.subplot(projection=wcscut, slices=('x','y',0,0))
    # im = ax.imshow(data)
    ax.contour(data, levels=[0.5, 1.0, 1.5, 2.0, 3.0], colors='red')
    # cbar = fig.colorbar(im)
    # cbar.set_label(r'Flux Density [Jy beam$^{-1}$]')

    # Display Herschel footprint
    spaxel = 9.4 / 3600.0
    # deciding what type of grayscale to make and then pulling out the values of the spaxels accordingly
    # by calling respective functions.
    if grayscale == "Abundance":
        # calls the abundance calculator code.
        import abundance_calculator as ac
        circle_intensities = ac.source_abundance(file)
        #if you want to plot the error gray scale you will need to use the following
        #circle_intensities = ac.source_abundance_error(file)
    elif grayscale == "Density":
        # calls the electron density calculator code.
        circle_intensities, _ = electrondensity.new_density(source) #[2592.049696547341, 3054.862181240571, 4069.1601640310564, 5613.428122094341, 12572.45617353372, 4828.190433699186, 5656.769219361611, 5707.226354181071, 5555.830259154675, 7903.5606254477325, 5728.5155269945, 5114.254914425617, 9058.242367674895, 3981.6723231222954, 3306.7077928167673, 10334.006436418536, 10293.870500235704, 6670.584456806155, 4651.879771630991, 2155.620648561293, 12845.536560461256, 12747.527582562903, 10570.833432168907, 4138.040923563868, 2478.6492178827993], [1041.3721796186212, 1282.8088385866354, 1927.3968340242866, 2787.495349013863, 5274.2475246335625, 2519.035163258442, 2867.8293309475143, 3040.2527323222903, 2442.0082962681436, 3617.783184568553, 3158.68396626751, 2722.972140753578, 3060.3313495092, 1839.3259990754632, 1456.492749055028, 5204.088165972358, 4222.853220184401, 3417.3854366415317, 2404.042012459343, 827.4970064963527, 6059.356606013506, 7750.490820369077, 5342.086313364377, 2142.084005209539, 987.5253463966739]
       # If you want to plot electron density error use the following.
        #_, circle_intensities = electrondensity.new_density(source)
        if transition_line != 88:
            return
    else:
        # calls the observed flux calculator code.
        circle_intensities = effecient.not_corrected_flux(file)

    declination_fit = f"/users/sdey/DATA/te_herschel/{file}_dec.fits"

    cube = fits.open(declination_fit)
    cube_header = cube[0].header
    cube_data = cube[1].data
    dec = cube_data[0]
    # print(dec)
    dec[[0, 4]] = dec[[4, 0]]
    dec[[1, 3]] = dec[[3, 1]]
    decData = dec.flatten()
    # dec_j2000 = np.round(decData, 4)
    # print("dec_j2000", dec_j2000)

    ra_fit = f"/users/sdey/DATA/te_herschel/{file}_ra.fits"

    cube = fits.open(ra_fit)
    cube_header = cube[0].header
    cube_data = cube[1].data
    ra = cube_data[0]
    ra[[0, 4]] = ra[[4, 0]]
    ra[[1, 3]] = ra[[3, 1]]
    raData = ra.flatten()

    increment = 0
    # list of pixel values which tells us the footprints pixel position on the radio image.
    x_val = []
    y_val = []

    cmap = cm.get_cmap('gray_r')
    norm = matplotlib.colors.Normalize(vmin=min(circle_intensities), vmax=max(circle_intensities))
    # calls the list of bad spaxels that will be skipped while plotting.
    if grayscale == "Density":
        skip_spaxels = electrondensity.bad_density_spaxels(source)
    else:
        skip_spaxels = sf.bad_spaxels(file)

    # Loops to draw the circles and numbers on the radio image. Also draws the grayscale.
    for i, (intensities, ra_val, dec_val) in enumerate(zip(circle_intensities, raData, decData)):
        x, y = wcs.world_to_pixel_values(raData[i], decData[i])
        x_val.append(x)
        y_val.append(y)
        if tuple(str(4 - increment) + str(i % 5)) in skip_spaxels:
            if i % 5 == 4:
                increment += 1
            continue
        r = SphericalCircle((ra_val * u.deg, dec_val * u.deg), spaxel / 2.0 * u.degree, edgecolor='black',
                            transform=ax.get_transform('fk4'), facecolor=cmap(norm(intensities)))
        ax.add_patch(r)
        ax.text(x, y, str(4 - increment) + ", " + str(i % 5), color="yellow", ha='center', va='center')
        if i % 5 == 4:
            increment += 1

    # making the color-bar.
    cbar_min = min(circle_intensities)
    # print(cbar_min)
    cbar_max = max(circle_intensities)
    # print(cbar_max)
    norm = plt.Normalize(vmin=cbar_max, vmax=cbar_min)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm)

    # shows colorbar units depending on the grayscale type.
    if grayscale == "Abundance":
        if transition_line == 88 or transition_line == 52:
            cbar.set_label(r'$\log(\frac{O^{++}}{H^+}) + 12$')
        elif transition_line == 57:
            cbar.set_label(r'$\log(\frac{N^{++}}{H^+}) + 12$')
        elif transition_line == 122:
            cbar.set_label(r'$\log(\frac{N^{+}}{H^+}) + 12$')
    elif grayscale == "Flux":
        cbar.set_label(r'FIR Flux (ergs s$^{-1}$ cm$^{-2}$)')
    else:
        cbar.set_label(r'Electron Density (cm$^{-3}$)')
    sm.set_cmap(cmap='gray_r')

    # save image
    # fig.tight_layout()
    plt.xlabel("RA (J2000)")
    plt.ylabel("Dec (J2000)")

    # print(source)
    if grayscale == "Density":
        if transition_line == 88:
            plt.title(f"{source}")
            plt.grid(True)
            plt.gca().set_aspect('equal')
            # print(min(raData), max(raData))
            # print(min(decData), max(decData))
            ax.set_xlim(min(x_val) - 15, max(x_val) + 15)
            ax.set_ylim(min(y_val) - 15, max(y_val) + 15)
            # if grayscale == "Flux":
            # plt.savefig(f"/users/sdey/DATA/flux_grayscale/{file}_flux")
            return plt.show()
    else:
        if transition_line == 88:
            plt.title(f"{source} [OIII] {str(transition_line)} microns ")

        elif transition_line == 52:
            plt.title(f"{source} [OIII] {str(transition_line)} microns ")

        elif transition_line == 57:
            plt.title(f"{source} [NIII] {str(transition_line)} microns ")
        else:
            plt.title(f"{source} [NII] {str(transition_line)} microns ")

        plt.grid(True)
        plt.gca().set_aspect('equal')
        # print(min(raData), max(raData))
        # print(min(decData), max(decData))
        ax.set_xlim(min(x_val) - 15, max(x_val) + 15)
        ax.set_ylim(min(y_val) - 15, max(y_val) + 15)
        #if grayscale == "Flux":
            #plt.savefig(f"/users/sdey/DATA/flux_grayscale/{file}_flux")
        return plt.show()

#not useful at the moment.
def cross_abundance_grayscale(file, Type):
    import abundance_calculator as ac

    file = file
    grayscale = Type
    source = file.split("_")[0]
    transition_line = int(file.split("_")[1])
    if transition_line != 88:
        return
    hdu = fits.open(f"/users/sdey/DATA/{source}_radio.fits")
    header = hdu[0].header
    data = hdu[0].data
    data = data[0, 0]
    wcs = WCS(header)
    # drop freq and pol axes
    wcs = wcs.dropaxis(3)
    wcs = wcs.dropaxis(2)

    ra_b1950, dec_b1950 = wcs.wcs.crval
    # Convert the B1950 coordinates to J2000
    coord_b1950 = SkyCoord(ra=ra_b1950 * u.deg, dec=dec_b1950 * u.deg, frame='fk4', equinox='B1950.0')
    coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))

    # Update the WCS to J2000
    wcs.wcs.crval = coord_j2000.ra.deg, coord_j2000.dec.deg
    ax = plt.subplot(projection=wcs)
    # ax = plt.subplot(projection=wcscut, slices=('x','y',0,0))
    # im = ax.imshow(data)
    ax.contour(data, levels=[0.5, 1.0, 1.5, 2.0, 3.0], colors='red')
    # cbar = fig.colorbar(im)
    # cbar.set_label(r'Flux Density [Jy beam$^{-1}$]')

    # Display Herschel footprint
    spaxel = 9.4 / 3600.0


    declination_fit = f"/users/sdey/DATA/te_herschel/{file}_dec.fits"

    cube = fits.open(declination_fit)
    cube_header = cube[0].header
    cube_data = cube[1].data
    dec = cube_data[0]
    # print(dec)
    dec[[0, 4]] = dec[[4, 0]]
    dec[[1, 3]] = dec[[3, 1]]
    decData = dec.flatten()
    # dec_j2000 = np.round(decData, 4)
    # print("dec_j2000", dec_j2000)

    ra_fit = f"/users/sdey/DATA/te_herschel/{file}_ra.fits"

    cube = fits.open(ra_fit)
    cube_header = cube[0].header
    cube_data = cube[1].data
    ra = cube_data[0]
    # print(ra)
    ra[[0, 4]] = ra[[4, 0]]
    ra[[1, 3]] = ra[[3, 1]]
    raData = ra.flatten()
    # ra_j2000 = np.round(raData, 4)
    # print("ra_j2000", ra_j2000)

    # coord_j2000 = SkyCoord(ra=ra_j2000 * u.deg, dec=dec_j2000 * u.deg, frame="fk5", equinox="J2000")
    # coord_b1950 = coord_j2000.transform_to(FK4(equinox="B1950"))
    # ra = coord_b1950.ra.deg
    # dec = coord_b1950.dec.deg
    # print(ra)
    # print(dec)
    increment = 0
    # list of pixel values which tells us the footprints pixel position on the radio image.
    x_val = []
    y_val = []
    # deciding what type of grayscale to make and then pulling out the values of the spaxels accordingly
    # by calling respective functions.
    if grayscale == "NN":
        # calls the abundance calculator code.
        circle_intensities, skip_spaxels = ac.nitrogen_ratio(source)
    else:
        circle_intensities, skip_spaxels = ac.nitroxy_ratio(source)
    cmap = cm.get_cmap('gray_r')
    norm = matplotlib.colors.Normalize(vmin=min(circle_intensities), vmax=max(circle_intensities))
    # calls the list of bad spaxels that will be skipped while plotting.


    # Loops to draw the circles and numbers on the radio image. Also draws the grayscale.
    for i, (intensities, ra_val, dec_val) in enumerate(zip(circle_intensities, raData, decData)):
        x, y = wcs.world_to_pixel_values(raData[i], decData[i])
        x_val.append(x)
        y_val.append(y)
        if tuple(str(4 - increment) + str(i % 5)) in skip_spaxels:
            if i % 5 == 4:
                increment += 1
            continue
        r = SphericalCircle((ra_val * u.deg, dec_val * u.deg), spaxel / 2.0 * u.degree, edgecolor='black',
                            transform=ax.get_transform('fk4'), facecolor=cmap(norm(intensities)))
        ax.add_patch(r)
        ax.text(x, y, str(4 - increment) + ", " + str(i % 5), color="yellow", ha='center', va='center')
        if i % 5 == 4:
            increment += 1

    # making the color-bar.
    cbar_min = min(circle_intensities)
    # print(cbar_min)
    cbar_max = max(circle_intensities)
    # print(cbar_max)
    norm = plt.Normalize(vmin=cbar_max, vmax=cbar_min)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm)

    # shows colorbar units depending on the grayscale type.
    if grayscale == "NN":
        cbar.set_label(r'$\log(\frac{N^{++}}{N^+})$')
    else:
        cbar.set_label(r'$\log(\frac{N^{++}}{O^{++}})$')

    sm.set_cmap(cmap='gray_r')
    # save image
    # fig.tight_layout()
    plt.xlabel("RA (J2000)")
    plt.ylabel("Dec (J2000)")

    # print(source)
    plt.title(f"{source}")
    plt.grid(True)
    plt.gca().set_aspect('equal')
    # print(min(raData), max(raData))
    # print(min(decData), max(decData))
    ax.set_xlim(min(x_val) - 15, max(x_val) + 15)
    ax.set_ylim(min(y_val) - 15, max(y_val) + 15)
    # if grayscale == "Flux":
    # plt.savefig(f"/users/sdey/DATA/flux_grayscale/{file}_flux")
    return plt.show()

