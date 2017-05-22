from Headers import*
      
class PythonPeaks:
    
    def __init__(self, obsid):
        
        self.obsid = obsid
              
        #Parameter used in peaks masking        
        self.width = 5 # [GHz]
        
        # Noise estimation
        self.sigma = 5    
        
    ''' Get level-2 spectrum form the folder.'''
    def getSpireFtsLevel2(self, spectrumType, tmpDir, what='spss'): # 'WHAT = 'sds' etended calibration
    
        tarFile = tmpDir + "/%i_FTS_sparse_level2.tar"%self.obsid
    
        with tarfile.open(tarFile,'r') as tar:
            for member in tar.getmembers():
                if (what in member.name and spectrumType in member.name):
                    f=tar.extractfile(member)
                    xx = fits.open(gzip.open(io.BytesIO(f.read())))
        spec = {}
        with xx as hdu:
            #
            for k in hdu:
                extname = k.name
                if ('S' in extname):
                    spec[k.name] = {}
                    spec[k.name]['wave'] = k.data["wave"]
                    spec[k.name]['flux'] = k.data["flux"]
                    spec[k.name]['fluxErr'] = k.data["error"]
     
        return spec
            
    ''' Fitting continuum based on coefficients form Ros' feature finder '''
    def fitContinuum(self, pixel, frequency, continuum):

        param = []          
        param.append(continuum[pixel]['p3'])
        param.append(continuum[pixel]['p2'])   
        param.append(continuum[pixel]['p1'])
        param.append(continuum[pixel]['p0'])
        continuumCoef = []
        degree = len(param) - 1
    
        for i in range (0, len(frequency)):
            y= param[degree]
            for j in range(0, degree):
                y += param[j] * frequency[i] ** (degree- j)
            continuumCoef.append(y)
            
        return continuumCoef
        
    ''' Substracting continuum from the spectra
         - returns array of substracted data '''
    def substractContinuum(self, flux, continuum):
            
        substracted= [] 
        for i in range (0, len(flux)):
            substracted.append(flux[i] - continuum[i])  
            
        return substracted
    
    ''' Noise estimation using wavelet function to discard biggest peaks '''    
    def waveletPeaksSearching(self, flux, freq):
        
        noise = 25 # [%]
        SNR = 8
        
        noiseFreq = []
        
        wt_scales = np.array([3,4,5,6,7,8],dtype=np.int64)
        peak_idx = signal.find_peaks_cwt(flux, wt_scales, wavelet=None, min_snr= SNR, noise_perc = noise) 
        
        for i in peak_idx: noiseFreq.append(freq[i])
            
#        fig = plt.figure(figsize=(9,6),dpi=300)  
#        ax = fig.add_subplot(1,1,1) 
#        ax.grid(True)
#        ax.plot(freq, flux, 'k-')
#        ax.plot(freq[peak_idx],flux[peak_idx],'ro')
#        ax.set_title(r'Peaks detection using wavelet function.')
#        ax.set_xlabel('Frequency [GHz]')
#        ax.set_ylabel('Flux [Jy]')
#        
        return noiseFreq
        
        
    ''' Flagging regions nearby the peaks +/- N GHz 
        - returns frequency, and flux arrays with masked regions'''
    def maskRegions(self, frequency, flux, peaks):
                                    
        newFreq = []       
        for i in frequency:
            for j in peaks:
                if np.abs(i - j) <= self.width:
                    newFreq.append(i)
     
        newFreq = np.unique(newFreq)
        
        for j in newFreq:             
            index = np.argwhere(frequency == j)
            frequency = np.delete(frequency, index) 
            flux = np.delete(flux, index) 
        
        noise =  {"frequency":frequency, "flux":flux}

        return noise
        

    """ Noise estimation
        - median absolute deviation from data with biggest peaks ommited """
    def noiseEstimationMedian(self, flux):
        
        noise = 1.4826 * median_absolute_deviation(flux)  
              
        return noise
        
    ''' Emission lines:
        - searching for the peaks above selected noise. Noise = 5 sigma in overlapping region, and 3 simga in the rest''' 
    def maxSearching(self, frequency, flux, noise, detector):
                
        peakFreq = []
        peakSNR = []
        peakFlux = []

        # Searching for local extremums
        localMaxIndex = argrelextrema(np.array(flux), np.greater)[0]

        # Discarding first local extremum for each detector
        if detector == 'SSW':
          skip = len(localMaxIndex) - 1
          detEdge = min(frequency)

        if detector == 'SLW':
          skip = 0
          detEdge = max(frequency)
               
        # Searching for peaks above the noise    
        #   - additional condition: detections +/- 5 GHz from the overlapping area are automatically discarded
        for i in localMaxIndex:                        
          if flux[i] >= noise[i] and flux[i] != flux[localMaxIndex[skip]] and np.abs(frequency[i] - detEdge) > 5:
            peakSNR.append(flux[i]/noise[i])
            peakFreq.append(frequency[i])
            peakFlux.append(flux[i])
            
        featuresFound =  {"frequency":peakFreq, "SNR":peakSNR, "flux":peakFlux}
     
        return featuresFound
            
    def getSpireFtsHpdp(self, obsid, hpdpDir=""):

        files = glob.glob("%s/%i*"%(hpdpDir,obsid))

        if (len(files) < 1):
            print ("No HPDP file found for OBSID %i"%obsid)
            return None
            
        hdu = fits.open(files[0])   
        spec = {}
        for k in hdu:
            extname = k.name
            if (('SLW' in extname) or ('SSW' in extname)):
                spec[k.name] = {}
                spec[k.name]['wave'] = k.data["wave"]
                spec[k.name]['flux'] = k.data["flux"]
                spec[k.name]['fluxErr'] = k.data["error"]
        return spec
    

     