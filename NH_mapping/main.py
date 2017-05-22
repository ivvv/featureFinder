from PP_classSpectralCubes import PythonPeaks
from Headers import *
from specCubes_NH import getSpireFtsCube

version = "1.0_IV" # 22/05/2017, IV modifications
# Continuum dir
tmpDir = home + "/featureFinder/mapping/HRmap/CP/continuum/"
#tmpDir = "D:\HRmap/CP/continuum/"
allFiles = os.listdir(tmpDir)

# Automatically create dir to save results. Creates "SpectralCubes folder in "home + Desktop/Dropbox/Work"
dir1 = home + "/featureFinder/mapping/PP_results/"
#dir1 = home + "/Desktop/Dropbox/Work/SpectralCubes/"
saveFits = dir1 + "PPspectralCubesFits"
saveCSV = dir1 + "PPspectralCubesCSV"
savePlots = dir1 + "PPspectralCubesPlots"

if not os.path.exists(dir1):
    os.makedirs(dir1)
    
if not os.path.exists(saveFits):
    os.makedirs(saveFits)
    
if not os.path.exists(saveCSV):
    os.makedirs(saveCSV)
    
if not os.path.exists(savePlots):
    os.makedirs(savePlots)

#%%
# Few obsids to try on
#obsid = 1342192174
#obsid = 1342245107
#obsid = 1342214846
#obsid = 1342214827

# Iterate over all observations
for i in range(0, len(allFiles)):
    
    start_time = time.time()
    obsid = int(allFiles[i][:10])

    # FITS file to save
    toFitsFreq = []
    toFitsError = []
    toFitsSNR = []
    toFitsFlux = []
    toFitsDetector = []
    toFitsRow = []
    toFitsColumn = []
    toFitsFlag = []
    
    #CSV file to save
    csv_file = open( saveCSV +'/%i_featuresFoundPP_%s.csv'%(obsid,version), 'w', newline = '')
    featuresFoundNH = csv.writer(csv_file)        
    featuresFoundNH.writerows([['frequency','frequencyError','SNR','fluxSub','detector', 'row', 'column', 'featureFlag']])
    featuresFoundNH.writerows([['Double','Double','Double','Double','String', 'String','String','Double']])
    featuresFoundNH.writerows([['GHz','GHz','None','Jy','None','None','None','None']])
    featuresFoundNH.writerows([['Measured frequency','Error on frequency','Peak/local noise','Subtracted continuum flux','Detector', 'Row', 'Column', 'Feature flag']])
        
    # Check if continuum for sepcific obsid exists
    try:
        csvFile = tmpDir + "%s_fittedContinuumParameters.csv"%obsid
        continuumFile = ascii.read(csvFile, data_start=4, comment = '#') 
    except:
        print('No continuum!')
        continue
   
    # Try to get apodised data from the spec cube    
    try:
        result = getSpireFtsCube(obsid, apod=True)
       
        if result == None:
            raise 
    except:
        print('No data to import!')
        continue 
#%%        
    # Detectors
    detector = ["SLW", "SSW"]
    # The target name
    target = result["SSW_HR_convol"].header["OBJECT"]
    
    print("Doing {} ({}/{}), target = {}".format(obsid,i+1,len(allFiles),target))
    #print(obsid)
    
    # Iterate over detectors
    for detectorName in detector:
        
        if detectorName == "SLW":
            colour = 'y'
        else:
            colour = 'b'
            
        # Shape
        shape = result["%s_HR_convol"%detectorName].shape
        # central pixel
        xc, yc = int(shape[1]/2),int(shape[2]/2)
        # Get pixel
        for x_pixel in np.arange(0, shape[1]):
            for y_pixel in np.arange(0, shape[2]):
               
                    #print('Pixel = (%i, %i)'%(x_pixel, y_pixel)) 
                    
                    cspec = {}
                
                    for array in [detectorName]:
                        
                        # Extracting pixel's flux and frequency
                        cspec[array] = result["%s_HR_convol"%array][:,x_pixel,y_pixel]
                        freq = (cspec[array].spectral_axis/1.0e9).value
                        flux = cspec[array].value
                  
                    # Try if and pixel is empty/if yes ommit it
                    try:
                        if np.isnan(flux).any() == True:
                            raise
                        if np.isnan(freq).any() == True:
                            raise   
                    except:   
                        #print('No data')
                        continue
                       
                    # Get index of the pixel (x,y) from the continuum file
                    try:
                       pixel = np.where((continuumFile['row'] == x_pixel) & (continuumFile['column'] == y_pixel) & (continuumFile['array'] == detectorName))[0][0]
                    except:
                        print('No continuum')
                    
                    peaks = PythonPeaks(obsid)
                    
                    # Fitting continuum from the fits files
                    continuum = peaks.fitContinuum(pixel, freq, continuumFile)
                   
                    # Substract spectra
                    subContinuum = peaks.substractContinuum(flux, continuum)
                    
                    # Find biggest peaks with the wavelet function
                    peaksFound = peaks.waveletPeaksSearching(flux, freq)
                  
                    # Discarding peaks (emission lines)
                    maskedRegion = peaks.maskRegions(freq, subContinuum, peaksFound)
                                      
                    # Noise estimation
                    madNoise = peaks.noiseEstimationMedian(maskedRegion['flux'])
                    
                    # Plotting noise for emission lines
                    #   - use sigma = 5 for the overlapping area (above 900 for SLW, and below 1000 for SSW), and sigma = 3 for the rest

                    if detectorName ==  "SLW":                       
                        noiseMain = np.array([3 * madNoise] * len(freq[np.where(freq < 900)]))
                        noiseOverlapping = np.array([5 * madNoise] * len(freq[np.where(freq >= 900)]))
                        noise = np.concatenate((noiseMain, noiseOverlapping), axis = 0)  
                        
                    if detectorName == "SSW":
                        noiseMain = np.array([3 * madNoise] * len(freq[np.where(freq > 1000)]))
                        noiseOverlapping = np.array([5 * madNoise] * len(freq[np.where(freq <= 1000)]))                     
                        noise = np.concatenate((noiseOverlapping, noiseMain), axis = 0)  
                        
                    # Peaks searching   
                    featuresFound = peaks.maxSearching(freq, subContinuum, noise, detectorName )
                    
                    #Error
                    error = ((max(freq) - min(freq))/len(freq))/2
                    ##################PLOTS####################################
                    
                    ##########################################################
                    # Print oryginal spectrum
#                    fig = plt.figure(figsize=(13,5))
#                    ax = fig.add_subplot(1,1,1)
#                    ax.set_title('%i:%i OBSID: %i'%(x_pixel, y_pixel,obsid))
#                    plt.plot(freq,flux)
##                    plt.savefig(savePlots +'/%i_specCubeOryginal_%i_%i.png'%(obsid, x_pixel, y_pixel), format='png', dpi=200)
#                    ###########################################################
#                    
#                    ###########################################################
                    # Plotting (sepctra, continuum, subtracted spectra):
#                    fig0 = plt.figure(figsize=(12,8), dpi=200)
#                    ax0 = fig0.add_subplot(1,1,1)
#                    ax0.grid(True)
##                    ax0.set_title('%i:%i OBSID: %i, original apodized spectrum.'%(x_pixel, y_pixel,obsid))
#                    ax0.set_title('%i:%i OBSID: %i'%(x_pixel, y_pixel,obsid))
#                    ax0.set_xlabel('Frequency [GHz]')
#                    ax0.set_ylabel('Flux')                    
#                    #   1. Spectra
#                    ax0.plot(freq, flux, 'y')                 
#                    #   2. Continuum 
#                    ax0.plot(freq, continuum, 'k') 
#                    #   3. Substracted spectra 
#                    ax0.plot(freq, subContinuum, colour)                     
#                    plt.savefig(savePlots +'/%i_specCubeContinuum_%s_%i_%i.png'%(obsid, detectorName, x_pixel, y_pixel), format='png', dpi=200)
#                   ############################################################
#                   
#                    ############################################################
#                    # Plotting (subtracted continuum, masked region, noise):
#                    fig = plt.figure(figsize=(12,8), dpi=200)
#                    ax = fig.add_subplot(1,1,1)
#                    ax.grid(True) 
#                    ax.set_title('%i:%i OBSID: %i, subtracted continuum and masked region for noise estimation.'%(x_pixel, y_pixel,obsid))
#                    ax.set_xlabel('Frequency [GHz]')
#                    ax.set_ylabel('SNR')                    
#                    #   1. Subtracted continuum
#                    ax.plot(freq, subContinuum, colour) 
#                    #   2. Masked region
#                    ax.plot(maskedRegion['frequency'], maskedRegion['flux'], 'c')                            
#                    #   3. Noise (MAD)
#                    ax.plot(freq, noise,'g', linestyle= '--')
#                    ##########################################################
#
#                    ########################################################## 


# UNCOMMENT TO SAVE THE PLOTS FOR EACH PIXEL
                    if (x_pixel == xc and y_pixel == yc):
                        #Plotting (subtracted continuum, noise, PP results):
                        fig = plt.figure(figsize=(12,8), dpi=200)
                        ax2 = fig.add_subplot(1,1,1)
                        ax2.grid(True)                     
                        #   1. Subtracted continuum
                        ax2.plot(freq, subContinuum, colour,label='Continuum subtacted spectrum')                   
                        #   2. Median absolute value noise
                        ax2.plot(freq, noise,'g', linestyle= '--',label='SNR threshold')                                     
                        NHDetections = len(featuresFound['frequency']) 
                        ax2.plot(featuresFound['frequency'] , featuresFound['flux'], 'ro',label="Features")
                        ax2.set_title('Pixel:%i:%i, OBSID: %i, target: %s, PP = %i'%(x_pixel, y_pixel, obsid, target, NHDetections))
                        ax2.set_xlabel('Frequency [GHz]')
                        ax2.set_ylabel('Intensity [W/m2/Hz/sr]')
                        plt.legend(loc=2)
                        plt.savefig(savePlots +'/%i_specCube_%s_%i_%i_%s.png'%(obsid, detectorName, x_pixel, y_pixel,version), format='png', dpi=200)
#                    ############################################################  


                    # CSV file
                    for i in range(0, len(featuresFound['frequency'])):
                        results = []
                        results.append(featuresFound['frequency'][i])
                        results.append(error)
                        results.append(featuresFound['SNR'][i])
                        results.append(featuresFound['flux'][i])
                        results.append(detectorName)
                        results.append(x_pixel)  
                        results.append(y_pixel)   
                        results.append(0)    
                        featuresFoundNH.writerows([results])   
                        
                        toFitsFreq.append(featuresFound['frequency'][i])
                        toFitsError.append(error)
                        toFitsSNR.append(featuresFound['SNR'][i])
                        toFitsFlux.append(featuresFound['flux'][i])
                        toFitsDetector.append(detectorName)
                        toFitsRow.append(x_pixel)
                        toFitsColumn.append(y_pixel)
                        toFitsFlag.append(0)

                    
    csv_file.close()
#%%
    # FITS file   
    col1 = fits.Column(name='frequency', format='D', unit = 'GHz', array=toFitsFreq)
    col2 = fits.Column(name='frequencyError ', format='D', unit = 'GHz', array=toFitsError)
    col3 = fits.Column(name='SNR', format='D', array=toFitsSNR)
    col4 = fits.Column(name='fluxSub', format='D', array=toFitsFlux)
    col5 = fits.Column(name='detector', format='4A', array=toFitsDetector)
    col6 = fits.Column(name='row', format='I', array=toFitsRow)
    col7 = fits.Column(name='column', format='I', array=toFitsColumn)
    col8 = fits.Column(name='featureFlag', format='D', array=toFitsFlag)
    
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(saveFits + '/%i_featuresFoundPP_%s.fits'%(obsid,version),overwrite=True)
    
    # Time for this obsid
    elapsed_time = time.time() - start_time
    print('Time = %.3f s'%elapsed_time)
    #break
print ("All done")
#
