#!/usr/bin/env python3
import argparse, sys, numpy, os
import matplotlib.pyplot
import datetimelib
import spectrumlib
import generallib
from astropy.io import fits



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Loads all of the JSON format spectra in a JSON file.')
    parser.add_argument('list', type=str, help='List of all of the spectra (filenames).')
    parser.add_argument('--noplot', action="store_true", help="Suppress plotting of the input spectra.")	
    parser.add_argument('-c', '--config', type=str, help="Configuration file for the trail plots. Default is 'trail.cfg'.")
    parser.add_argument('-e', '--ephem', type=str, help="Ephemeris file.")
    parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
    parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
    parser.add_argument('-b', '--boost', action="store_true", help="Boost the trail contrast.")
    arg = parser.parse_args()
    hasEphem = False
    
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
        spectra.append(spectrum)
        #print("%d : %s"%(index+1, spectrum.name))
    print("Loaded %d spectra."%len(spectra))

    for s in spectra:
        s.phase = ephem.getPhase(s.HJD)
        print(s.HJD, s.phase)

    spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
	
    for s in spectra:
        print(s.HJD, s.phase)


    lowerWavelength = 4500
    upperWavelength = 6660
    print("Trimming to %f - %f Angstrom"%(lowerWavelength, upperWavelength))
    # Trim to the required wavelength range
    for s in spectra:
        #print(s.wavelengthRange)
        #print(len(s.wavelengths), len(s.flux), len(s.fluxErrors))
        if len(s.fluxErrors)==0:
            s.fluxErrors = numpy.ones(len(s.wavelengths))
        s.trimWavelengthRange(lowerWavelength, upperWavelength)

	# Resample down to lowest resolution
    minElements = 1E6
    for index, s in enumerate(spectra):
        elements = len(s.flux)
        if elements < minElements:
            minElements = elements
            shortest = index
	
    print("Shortest spectrum is number %d with %d elements."%(shortest, minElements))
    sampleWavelengths = spectra[shortest].wavelengths
    print(sampleWavelengths)
    print("Resampling...")
    for s in spectra:
	    s.resample(sampleWavelengths)

    
    totalAdded = 0
    numPhaseBins = 30
    ySize = numPhaseBins * 2
    xSizeArray = [len(s.wavelengths) for s in spectra]
    xSize = max(xSizeArray)
    trailArray = []
    phasedArray = []
    
    for index in range(ySize):
        blankSpectrum = numpy.zeros(xSize)
        trailArray.append(blankSpectrum)
        phasedArray.append(blankSpectrum)

    phaseLabels = []
    for bin in range(ySize):
        phaseLower = 2.0/ySize * bin
        phaseUpper = 2.0/ySize * (bin+1)
        phaseLabels.append((phaseLower+phaseUpper)/2)
        spectraToAdd = []
        for s in spectra:
            if s.phase<phaseUpper and s.phase>=phaseLower: 
                spectraToAdd.append(s)
            if (s.phase+1)<phaseUpper and (s.phase+1)>=phaseLower: 
                spectraToAdd.append(s)
        totalSpectrum = numpy.zeros(xSize)
        for s in spectraToAdd:
            totalSpectrum+=s.flux
        averageSpectrum = totalSpectrum / len(spectraToAdd)
        phasedArray[bin] = averageSpectrum
        print("Bin #%d  phase range %f to %f contains %d spectra to average together."%(bin, phaseLower, phaseUpper, len(spectraToAdd)))
        totalAdded+= len(spectraToAdd)

        print("%d spectra added to %d bins"%(totalAdded, ySize))


		

    # Set up the matplotlib environment
    generallib.setMatplotlibDefaults()
    params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
    matplotlib.rcParams.update(params)

    xValues = spectra[0].wavelengths
    for i, p in enumerate(phasedArray):
        yValues = p
        print("phase: %.2f"%phaseLabels[i])
        matplotlib.pyplot.xlim(lowerWavelength, upperWavelength)
        matplotlib.pyplot.ylim(-5, 20)
        matplotlib.pyplot.step(xValues, yValues)
        matplotlib.pyplot.text(5000, 18, "%.2f"%phaseLabels[i])
        matplotlib.pyplot.draw()
        matplotlib.pyplot.show(block=False)
        matplotlib.pyplot.pause(0.1)
        matplotlib.pyplot.clf()
    

    sys.exit()

