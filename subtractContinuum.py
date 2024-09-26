#!/usr/bin/env python3
import argparse, sys, numpy, os, json, copy
import matplotlib.pyplot
import datetimelib
import spectrumlib
import generallib
from astropy.io import fits


class configClass():
    def __init__(self):
        self.filename = ""
	
    def load(self, filename="flatten.json"):
        try:
            fileHandle = open(filename, 'rt')
            self.filename = filename
            configObject = json.load(fileHandle)
        except FileNotFoundError:
            print("Could not open config file %s"%filename)
            return
		
        for key in configObject.keys():
            keyString = str(key)
            value = configObject[key]
            if isinstance(value, (str)): 
                value = str(value)
            if isinstance(value, (list)):
                value = numpy.array(value)
            setattr(self, key, value)

        fileHandle.close()
		



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Loads all of the JSON format spectra in a JSON file.')
    parser.add_argument('list', type=str, help='List of all of the spectra (filenames).')
    parser.add_argument('--noplot', action="store_true", help="Suppress plotting of the input spectra.")	
    parser.add_argument('-c', '--config', default="flatten.json", type=str, help="Configuration file for the trail plots.")
    parser.add_argument('-e', '--ephem', type=str, help="Ephemeris file.")
    parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
    parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
    parser.add_argument('-b', '--boost', action="store_true", help="Boost the trail contrast.")
    arg = parser.parse_args()
    hasEphem = False
	
    config = configClass()
    config.load(arg.config)

    
    inputList = open(arg.list, "rt")
    fileList = []
    for line in inputList:
        line = line.strip()
        if line[0]=="#": continue
        if len(line)<4: continue
        fileList.append(line)


    if arg.ephem is not None:
        ephem = datetimelib.ephemerisObject()
        ephem.loadFromFile(arg.ephem)
        print("Loading ephemeris data: ", ephem)
        hasEphem = True
	
    spectra = []
    for index, filename in enumerate(fileList):
        spectrum = spectrumlib.spectrumObject()
        spectrum.loadFromJSON(filename)
        spectrum.loadedFrom = filename
        spectrum.name = os.path.splitext(filename)[0].split("/")[-1]
        if hasEphem: spectrum.phase = ephem.getPhase(spectrum.HJD)
        spectra.append(spectrum)
        #print("%d : %s"%(index+1, spectrum.name))
    print("Loaded %d spectra."%len(spectra))

    # Remove any negative values in the spectra
    for s in spectra:
        s.removeNegatives()
        s.removeZeros()
	
    
    # Set up the matplotlib environment
    generallib.setMatplotlibDefaults()
    params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
    matplotlib.rcParams.update(params)


    
    mask = [ [m['lower'], m['upper']] for m in config.mask ]
    print(mask)
    mask = spectrum.invertMask(mask)
    print(mask)

    # Create figure and axes
   
    # Show the mask for the continuum fits
    fig, ax = matplotlib.pyplot.subplots(figsize=(11, 11/1.618))
    for n, spectrum in enumerate(spectra):
        print(n+1, spectrum.loadedFrom)
        xValues = spectrum.wavelengths
        yValues = spectrum.flux

        matplotlib.pyplot.step(xValues, yValues)
        ax = matplotlib.pyplot.gca()
        yLims = ax.get_ylim()
        height = yLims[1] - yLims[0]
        for maskArea in mask:
            width = maskArea[1]-maskArea[0]
            rect = matplotlib.patches.Rectangle((maskArea[0], yLims[0]), width, height, facecolor='red', alpha=0.4)
            ax.add_patch(rect)
        
        yValues = spectrum.fitPoly(degree= 9, mask = mask)

        continuumSubtracted = spectrum.flux - yValues

        # Create a new spectrum and save it
        cSpectrum = copy.deepcopy(spectrum)
        cSpectrum.name = "continuum subtracted " + spectrum.name
        cSpectrum.flux = continuumSubtracted
        print(cSpectrum.name)

        cSpectrum.writeToJSON("continuum/" + cSpectrum.name)
        matplotlib.pyplot.plot(xValues, yValues)
        matplotlib.pyplot.plot(xValues, continuumSubtracted)
        
        matplotlib.pyplot.show(block = False)
        matplotlib.pyplot.pause(0.001)
        matplotlib.pyplot.clf()
        

    sys.exit()
    
    
    
    
    
    
    
    