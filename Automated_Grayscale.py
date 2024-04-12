import GrayScale_2000
import source_file

for name in ["W3_88","W3_122", "W3_52", "W3_57"]: #source_file.source_name():
    """Automates the production of grayscale"""
    try:
        """Calls the Grayscale function. The second parameter changes according to the grayscale type: "Density", "Flux", "Abundance"""
        GrayScale_2000.grayscale(name, "Density")
        #for i in ["NN","NO"]:
            #GrayScale_2000.cross_abundance_grayscale(name,i)
    except FileNotFoundError:
        continue

