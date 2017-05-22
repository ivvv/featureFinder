import os

################################################################################
####        FTS Spectral Feature Finder v22Final mapping PG 19/11/16        ####
################################################################################

# Script to find and fit features in SPIRE FTS spectra
#
# The script is set out in the following sections
###### User options ######
# 1. Example obsids for testing
# 2. Input options
# 3. Gof/ToFE options
# 4. Product save options
# 5. Plotting options
######## Functions #######
##Compiled from a separate file
# 6. Import statements
# 7. Plotting
# 8. Products I/O
# 9.  GoF/ToFE
# 10. Check fit/redshift
# 11. Feature finding
### Running the script ###
# 12. script execution

################################################################################
####                    1. Example obsids for testing                       ####
################################################################################

###########################
##Orion bar (full)
#obsids = [0x5000C8FE]

##IRAS18507+0121-1 (intermediate)
#obsids = [0x50016A5C]

##Quick intermediate
#obsids = [1342202260,1342202260] 
###########################

######Get HR mapping obsids from the obs log
obsLog = SpireObsLogProduct()
#SpireSpectroPoint + Raster mode + HR/CR/H+LR + intermediate + full
mode = obsLog["obsLog"]["obsMode"].data
res = obsLog["obsLog"]["specRes"].data
sampling = obsLog["obsLog"]["specSampl"].data
obsid = obsLog["obsLog"]["obsId"].data
modePos = mode.where((mode == "SpireSpectroPoint").or(mode == "Raster").and(res != "L").and(res != "M").and(sampling != "sparse"))
obs = obsLog["obsLog"].select(modePos)
odNum = obs["odNum"].data
ix = odNum.where(odNum >= 209)
nx = ix.length()
print "Running the feature finder on %i HR mapping/raster observations with OD >= 209"%nx
print "Although skipping a few (dark sky, Uranus)" #, OMC-1)"

obsSel = obs.select(ix)
#
skipTargets = (["DARKSKY","URANUS"]) #,"OMC-1","SGRA_F1","SGRA_F2","SGRA_F3","SGRA_F4","SGRA_F5","SGRA_F6","SGRA_F7"])
#
targets = obsSel["source"].data
obsids = obsSel["obsId"].data
targets = ['OrionBar', 'otherOne']
obsids = [1342214846, 1342214827]
targets = ['otherOne']
obsids = [1342214827]
###These should be skipped and run separatly (takes a week)
obsidsToSkip = [1342192173,1342192174,1342192175,1342204898,1342204920,1342214827,1342214841,\
                1342214846,1342228703,1342243631,1342243632,1342243633,1342243634,1342243635,\
                1342243636,1342243637,1342245851,1342262908,1342262909,1342262913,1342262916,\
                1342265845]
obsidsToSkip = []              
#
# GRID
#

nx = len(obsids)
if (Configuration.hasProperty('start')):
   startIdx = int(Configuration.getProperty('start'))
else:
   startIdx = 0
#
# end index
#
if (Configuration.hasProperty('end')):
   endIdx = int(Configuration.getProperty('end'))
else:
   endIdx = startIdx + 300
if (endIdx > nx): endIdx = nx
#startIdx = 600
#endIdx = nx
idx0 = startIdx
idxn = endIdx
#
# GRID
#

onGrid = True

if onGrid:
    runFolder = 'hrMappingPPtest'
else:
    runFolder = 'hrSparseV24Public_mappingTest'

################################################################################
####                             2. Input options                           ####
################################################################################

#Set useHsa to fetch observations from the archive
#Or HIPE will look for a pool
useHsa = True

#CP or Naive cubes?
#cubeName = 'cube_convol'
cubeName = 'cube_convol'

if cubeName == 'cube_convol':
    folderName = 'Cp'
else:
    folderName = 'Naive'

#Set the resolution of the observation(s)
res = 'HR'

apod = False

#Set printToScreen if you want all features printed to the console
printToScreen = False
#Set if you want the number of features found per spaxel printed to the console
printCubeVerbose = False

#Set the overlap region
overlap = [944.0,1017.8]
#The additional noisy regions (GHz), based on 
#HR sensitivity > 0.25 Jy
#HR overlap > 0.45 Jy
noisyRegions = [608.57868974,1504.6583467]

###FINDING PARAMETERS###
snrCut = 5.0 #final cut using final noise estimate from full residual
snrCap = 500. #Highest SNR allowed for the initial guess
snrRange = [4.0,25.0] #(GHz) range either side of feature for final noise estimate
avoid = 10 #(GHz) excludes the ends of the bands by this width
mergeRadius = 20.0 #(GHz) merge radius for preliminary features found
#flagWidths = [6.0, 4.0, 2.0, 2.0, 2.0, 2.0] #(GHz) HALF the width of region flagged for feature found
#snrThresholds = [100.0, 50.0, 30.0, 10.0, 5.0, 3.0] #SNR thresholds for finding to loop over
flagWidths = [8.0, 8.0, 5.0, 4.0, 2.0, 2.0] #(GHz) HALF the width of region flagged for feature found
snrThresholds = [100.0, 50.0, 30.0, 10.0, 5.0, 1.0] #SNR thresholds for finding to loop over
signs = ['pos']*len(snrThresholds)
#signs+= ['neg']*len(snrThresholds)
#flagWidths+=flagWidths
#snrThresholds+=snrThresholds
if apod:
    flagWidths = [6.0, 6.0, 4.0, 4.0, 2.0, 2.0, 0.0] #(GHz) HALF the width of region flagged for feature found
    snrThresholds = [100.0, 50.0, 30.0, 10.0, 5.0, 3.0]#, 1.0] #SNR thresholds for finding to loop over
    snrThresholds = [x/3. for x in snrThresholds]
    snrCut = 4
#Default colours for plotting per snrThreshold are [black,lightOrange,darkGreen,purple,cyan,teal]
polyOrder = 3 #order of polynomial fitted to the continuum
#Initial continuum fit only parameters
resampleResolution = 5.0 #(GHz) for continuum fitting only
jumpThreshold = 3.5 #minimum jump (in RMS) to flag for continuum fitting only
baselineOrder = 3 #order of poly used for final SNR estimate
gOfFRng = 3 #range used by GofF to take the stats

limitValue = 2. #To prevent features racing off or jumping outside of the frequency bands
checkValue = 2. #To remove features moved too far

##SET THE MAXIMUM NUMBER OF ITERATIONS FOR THE SPECTRUM FITTER
#maxIter=None
maxIter = 500

##Set checkCutResults to a snrThreshold if you want to check the corresponding 
##mid-cut residual and totalModel *produced for 1st detector in 1st obsid only*
checkCutResults = False
#checkCutResults = snrThresholds[3]

###LR needs a different set of some finding parameters
#if res == 'LR':
#    resampleResolution = 1.0
#    snrThresholds = [50.0, 30.0, 10.0, 5.0]
#    flagWidths = [30.0, 25.0, 15.0, 15.0]
#    snrRange = [4.0, 50.0]
#    mergeRadius = 150.0
#    avoid = 30.0
#    baselineOrder = 1 #order of poly used for final SNR estimate
#    gOfFRng=10.0
#    polyOrder = 2
#    #Removes features in the avoid region and outside of the band
#    #Removes features that have traveled too far
#    #Limits the feature position during the fit
#    extraConstraints=True
#    useSensitivityCurve = False #this is not working well yet
#    limitValue = 20
#else:
useSensitivityCurve = False
extraConstraints=False

#Set to test improving the initial parameter guess
#Not set up for LR
#if res != 'LR':
##################
##1. Check the initial peak position
testCheckFit = 0
##2. Check both spectral max and min
testMaxMinAmp = 0
##3. Set to attempt an estimate of the source redshift
##Set to: 0 to turn off; 1 to use redshift1 function; 2 to use redshift2 function
estimateRedshift = 2
##4. & 5. Set to test notoriously blended lines with two sinc profiles
testNeutralC = 0
##################
#maxFitThresh = 30.0 #SNR above which the check will be performed
maxFitThresh = 5.0 #SNR above which the check will be performed
sampleRange = 3.0
window = 20.0
velocityIter = Double1d([6000.0])
ALLFreqs = Double1d([461.0407682,576.2679305,691.4730763,806.6518060,921.7997000,\
            1036.9123930,1151.9854520,1267.0144860,1381.9951050,1496.9229090])
#array with the difference patern for 12CO lines
refArray = ALLFreqs - ALLFreqs[0]
if testMaxMinAmp:
    testCheckFit = 1
if testNeutralC:
    testCheckFit = 1
    SLWFreqs = Double1d([461.0407682,576.2679305,691.4730763,806.6518060,921.7997000])
    SSWFreqs = Double1d([1036.9123930,1151.9854520,1267.0144860,1381.9951050,1496.9229090])
#else:
#    testCheckFit = 0
#    testMaxMinAmp = 0
#    testNeutralC = 0
#    estimateRedshift = 2

if onGrid:
    pathToFunctions = '/home/herspire/herspire/SPEC/ftsFeatureFinder/'
    path = '%s%s/'%(pathToFunctions,runFolder)
else:
    pathToFunctions = '/Users/ros/spec_out/ftsSpectralFeatureFinder/scriptsEditedForTwiki/postV11/esacVersions/nextRound/nextRound/2ndRelease/1stPublicRelease/'
    path = '/Users/ros/testing/%s/'%(runFolder)

#log events and parameters used
date = '%s'%(java.util.Date())
###logfile name
#logFileName = '/path/to/log/file/FSFF_logfile_%s%s%s_%sh%sm%ss.txt'%(date[8:10],date[4:7],date[24:],date[11:13],date[14:16],date[17:19])
logFileName = '%slogs/FSFF_logfile_%i_%i_%s%s%s_%sh%sm%ss.txt'%(path,idx0,idxn,date[8:10],date[4:7],date[24:],date[11:13],date[14:16],date[17:19])
###set to True to write the log to a text file called logFileName
writeLogToFile = True

################################################################################
####                             3. Gof/ToFE options                        ####
################################################################################

#Now turned on as default#

###GoF
testGofFfinal = 1
testGofFcut = 0
testGofFcutLines=0

###ToFE
#Turned off for mapping for now
#as most of the results are NaN
testFitEvidence = 0
limitFitEvidence = True
plotFitEvidence = False
verboseFitEvidence = False
useWeights = False

################################################################################
####                          4. Product save options                       ####
################################################################################

###FEATURES FOUND###

###set an output directory to write the feature tables to
outDir = '/where/the/tables/are/saved/'
outDir = path
###TABLE PER OBSID###
###Save a table per obsid in a csv file
writeResultsObsid = 1

###COMBINED TABLE FOR ALL OBSID###
###Save the results for all obsids in one csv file
writeResultsCombined = 1

###FITTED CONTINUUM###

###Set an output directory for the saved continuum/continua
###Two cubes saved per obsid
continuumPath = '/where/the/continua/are/saved/'
continuumPath = '%scontinuum'%(path)
###Save a table of fitted parameter per obsid
###Written to csv.
saveContinuumTable = 1
###Save the fitted continuum curve per obsid
###An SDS written to FITS
saveContinuumFits = 0

################################################################################
####                           5. Plotting options                          ####
################################################################################

##**ONLY POSTCARDS PRODUCED FOR MAPS**##

###POSTCARDS###

###Set plotPostCards to plot a postcard per obsid
plotPostCards = 1
###Set to save the postcard plots
savePostCards = 1
###Set to close the postcard plots
#(This is not linked to savePostCards)
closePostCards = 1
###Set a folder to save the postcards to
postCardPlotPath = '/where/postcard/plots/are/saved/'
#postCardPlotPath = '%spostcards/'%(path)
postCardPlotPath = '%spostcards'%(path)

###If closePostCards is True, then the postcard plots are not shown during the script execution
if closePostCards:
    seePostCards = 0
else:
    seePostCards = 1

################################################################################
####                             6-11 functions                             ####
################################################################################
if onGrid:
    execfile('%sscripts/featureFinderFunctionsV24Final_mappingTest_forMapping.py'%(pathToFunctions))
else:
    execfile('%sfeatureFinderFunctionsV24Public_mappingTest_forMapping.py'%(pathToFunctions))
#execfile('%sscripts/featureFinderFunctionsV22FinalFor1stRelease.py'%(pathToFunctions))
#execfile('/Users/ros/spec_out/ftsSpectralFeatureFinder/scriptsEditedForTwiki/postV11/esacVersions/nextRound/nextRound/2ndRelease/1stPublicRelease/featureFinderFunctionsV24Public_mappingTest_forMapping.py')

################################################################################
####                          12. script execution                          ####
################################################################################

if res == 'LR':
    threshColours = threshColours[0:len(snrThresholds)]
plotCol = {}
for go in range(len(snrThresholds)):
    plotCol[snrThresholds[go]]=threshColours[go]

####Setting up the log file with number of observations
####options chosen and parameters used
if writeLogToFile:
    file = open(logFileName,'w')
    file.write('-----------------\n')
    file.write('   %s %s obsids\n'%(len(obsids),res))
    file.write('-----------------\n')
    file.write('### Cube type selected ###\n')
    file.write('cubeName:%s\n'%(cubeName))
    file.write('-----------------\n')
    file.write('### Parameters used ###\n')
    file.write('snrCut:%s, snrCap:%s, snrRange:[%s,%s]\n'%(snrCut,snrCap,snrRange[0],snrRange[1]))
    file.write('avoid:%s, mergeRadius:%s, limitValue:%s, checkValue:%s \n'%(avoid,mergeRadius,limitValue,checkValue))
    if maxIter:
        file.write(', maxIter:%s\n'%(maxIter))
    else:
        file.write('\n')
    file.write('flagWidths:[%s'%(flagWidths[0]))
    for flag in range(1,len(flagWidths)): file.write(', %s'%flagWidths[flag])
    file.write(']\n')
    file.write('snrThresholds:[%s'%(snrThresholds[0]))
    for thresh in range(1,len(snrThresholds)): file.write(', %s'%snrThresholds[thresh])
    file.write(']\n')
    file.write('signs:[%s'%(signs[0]))
    for sign in range(1,len(signs)): file.write(', %s'%signs[sign])
    file.write(']\n')
    file.write('polyOrder:%s, resampleResolution:%s,jumpThreshold:%s,\n'%(polyOrder,resampleResolution,jumpThreshold))
    file.write('baselineOrder:%s\n'%(baselineOrder))
    file.write('-----------------\n')
    file.write('### Testing options chosen ###\n')
    file.write('testCheckFit:%s, testMaxMinAmp:%s, testNeutralC:%s, estimateRedshift:%s \n'\
         %(testCheckFit,testMaxMinAmp,testNeutralC,estimateRedshift)\
             +'testGofFfinal:%s, testGofFcut:%s, testFitEvidence:%s\n'\
                 %(testGofFfinal,testGofFcut,testFitEvidence))
    #if res == 'LR': file.write('useSensitivityCurve:%s\n'%(useSensitivityCurve))
    file.write('-----------------\n')
    file.write('### Saving options chosen ###\n')
    file.write('writeResultsObsid:%s, writeResultsCombined:%s\n'\
         %(writeResultsObsid, writeResultsCombined))
    file.write('saveContinuumTable:%s, saveContinuumFits:%s\n'%(saveContinuumTable, saveContinuumFits))
    file.write('-----------------\n')
    file.write('### Plotting options chosen ###\n')
    #file.write('plotIt:%s, saveMainPlots:%s, zoomPlots:%s\n'%(plotIt,saveMainPlots,zoomPlots))
    file.write('plotPostCards:%s, savePostCards:%s\n'%(plotPostCards,savePostCards))
    file.write('-----------------\n')
    file.write('### Features found ###\n')
    file.close()

#Speed up loading observations
archive = HsaReadPool()
hsa_storage = ProductStorage([archive])
#hsa_storage.register(HsaReadPool)

arrays = ['SLW','SSW']

#Empty lists to fill for the combine results for all obsids
freqDetsAll, indDetsAll ,freqErrDetsAll ,snrDetsAll ,detDetsAll = [],[],[],[],[]
obsidAll, raAll ,decAll ,resolutionAll ,biasModeAll ,mapSamplingAll = [],[],[],[],[],[]
rowAll, columnAll = [], []
zEstimateAllAve, zEstimateAllSlw, zEstimateAllSsw, finalZall, finalZerrAll, finalZflagAll, nMaxMatchAll = [],[],[],[],[],[], []
if testGofFfinal:
    gOfFAll, gOfFflagsAll = [], []
    combinedFlagsAll, combinedFlagsArrayAll = [], []
    resultAll = None
    
if testFitEvidence:
    oddsAll, resultAll, fitEvidenceFlagsAll = [], [], []

###Loop over the obsids
for go in range(startIdx,endIdx):
#for go in range(obsids.size):
    itarget = re.sub(r"\s+", '',targets[go]).upper()
    if (itarget in skipTargets):
        message = '%s: skipping %i %s'%(date[11:19],obsids[go], itarget)
        if writeLogToFile:
            addToLog(message,logFileName)
        continue
    obsid = obsids[go]
    if obsid in obsidsToSkip:
        message = '%s %s: skipping %i %s'%(date[11:19],go,obsid, itarget)
        if writeLogToFile:
            addToLog(message,logFileName)
        continue
    #check to see if the catalogue for this obsid already exists
#    try:
#        checkFile = fitsReader('%scataloguesFits/%s_featuresFound.fits'%(outDir,obsid))
#        continue
#    except:
#        pass
    #log events
    date = '%s'%(java.util.Date())
    message =  '%s %s: '%(date[11:19],go)
    #Empty lists to fill and combine results per obsid for the output product
    freqArrayPix, initialFreqArrayPix, indArrayPix, freqErrArrayPix, snrArrayPix = [],[],[],[],[]
    arrayArrayPix, rowArrayPix, columnArrayPix, raArrayPix, decArrayPix, modelsArrayPix = [],[],[],[],[],[]
    snrIterArrayPix, threshArrayPix, residArrayPix, sincArrayPix, totalModArrayPix = [],[],[],[],[]
    expectedFreqPix = []
    nFeatures = Long1d([0,0])
    zEstimateDets, featuresCentDets, errCentDets, snrCentDets = [], [], [], []
    zEstimateAve = 0.0
    if testGofFfinal:
        gOfFArrayAll, gOfFflagsArrayAll = [], []
        combinedFlagArrayAll = []
        resultArrayAll = None
    if testFitEvidence:
        oddsArrayAll, resultArrayAll, fitEvidenceFlagsArrayAll = [], [], []
    ###Load the observation context and retrive the necessary product
    obs, product, name = getProduct(obsid, res=res, useHsa=useHsa, mapping=True, cubeName=cubeName)
    #ppTableIn = fitsReader('/Users/ros/spec_out/ftsSpectralFeatureFinder/python/natalia/mappingTestMay17/firstTest/PPmapping/PPspectralCubesFits/%s_featuresFoundPP.fits'%(obsid))
    ppTableIn = fitsReader('/Users/ros/spec_out/ftsSpectralFeatureFinder/python/natalia/mappingTestMay17/firstTest/PPmapping/PPspectralCubesFits/%s_featuresFoundPP.fits'%(obsid))
    message+='%s, '%(name)
    if product[0] == 'failed':
        print "0x%X failed, skipping it"%(obsid)
        if writeLogToFile:
            addToLog(message,logFileName)
        continue
    print '\n%s 0x%X %s (index %i %i/%i)'%(obs.meta['object'].value,obsid,obsid,go,go-idx0+1,idxn-idx0)
    #create the table to fill with the fitted polynominal parameters
    continuumTable = TableDataset(description = 'Polynomial parameters fitted to the continuum, of the form p0 + p1*freq + p2*freq**2 + p3*freq**3')
    continuumTable.meta['obsid']=LongParameter(obsid)
    #create a structure to put the continuum parameters in during the array loop
    continuumCat = {}
    continuumCat['array']=String1d()
    continuumCat['row']=Int1d()
    continuumCat['column']=Int1d()
    continuumCat['ra']=Double1d()
    continuumCat['dec']=Double1d()
    continuumCat['param']=[]
    continuumCat['error']=[]
    #NEW FLAGGING
    flags = {'SLWC3':{'p1':Double1d(),'fL':Double1d(),'fR':Double1d()},\
             'SSWD4':{'p1':Double1d(),'fL':Double1d(),'fR':Double1d()}}
    cSpax = ' '
    for getSpecGo in range(len(product)):
        array = arrays[getSpecGo]
        cSpax+=array
        cube = product[getSpecGo]
        featuresArray = []
        name = '%s %s'%(cube.meta['object'].value,obsid)
        #make a copy of the input sds to fill with the final continuum fit.
        continuum = cube.copy()
        if plotPostCards:
            if array == 'SLW':
                plot2 = PlotXY(seePostCards) 
                #find the central pixel
                rowC = rowCslw = cube.wcs.naxis1/2
                colC = colCslw = cube.wcs.naxis2/2
                specSlwC = cube.getSpectrum1d(rowC, colC)
                rowC = rowCssw = product[1].wcs.naxis1/2
                colC = colCssw = product[1].wcs.naxis2/2
                specSswC = product[1].getSpectrum1d(rowC, colC)
                #check for NaNs
                if specSlwC.flux.where(IS_NAN(specSlwC.flux)).toInt1d().size > 0:
                    rowCslw = None
                    for rowC in range (cube.dimensions[1]):
                        for colC in range (cube.dimensions[2]): 
                            specSlwC = cube.getSpectrum1d(rowC, colC)
                            if specSlwC.flux.where(IS_NAN(specSlwC.flux)).toInt1d().size == 0:
                                rowCslw, colCslw = rowC, colC
                                break
                        if rowCslw != None:
                            break
                if specSswC.flux.where(IS_NAN(specSswC.flux)).toInt1d().size > 0:
                    rowCssw = None
                    for rowC in range (product[1].dimensions[1]):
                        for colC in range (product[1].dimensions[2]):
                            specSswC = product[1].getSpectrum1d(rowC, colC)
                            if specSswC.flux.where(IS_NAN(specSswC.flux)).toInt1d().size == 0:
                                rowCssw, colCssw = rowC, colC
                                break
                        if rowCssw != None:
                            break
                rowCs = [rowCslw,rowCssw]
                colCs = [colCslw,colCssw]
                #Create shaded regions
                boxMinMax = [MIN([MIN(specSlwC.flux),MIN(specSswC.flux)]),\
                              MAX([MAX(specSlwC.flux)*2,MAX(specSswC.flux)*2])]
                boxX = Double1d([0,0,avoid,avoid,0])
                boxY = Double1d([boxMinMax[0]-boxMinMax[1],boxMinMax[1],boxMinMax[1],boxMinMax[0]-boxMinMax[1],boxMinMax[0]-boxMinMax[1]])
                layB1 = LayerXY(boxX+specSlwC.wave[0],boxY,stroke=0,color=Color(0,0,0,0))
                layB1.style.fillClosureType=FillClosureType.SELF
                layB1.style.fillPaint=fillColourSlw
                layB1.style.fillEnabled=1
                layB2 = LayerXY(boxX+specSlwC.wave[-1]-avoid,boxY,stroke=0,color=Color(0,0,0,0))
                layB2.style.fillClosureType=FillClosureType.SELF
                layB2.style.fillPaint=fillColourSlw
                layB2.style.fillEnabled=1
                layB3 = LayerXY(boxX+specSswC.wave[0],boxY,stroke=0,color=Color(0,0,0,0))
                layB3.style.fillClosureType=FillClosureType.SELF
                layB3.style.fillPaint=fillColourSsw
                layB3.style.fillEnabled=1
                layB4 = LayerXY(boxX+specSswC.wave[-1]-avoid,boxY,stroke=0,color=Color(0,0,0,0))
                layB4.style.fillClosureType=FillClosureType.SELF
                layB4.style.fillPaint=fillColourSsw
                layB4.style.fillEnabled=1
                plot2.addLayer(layB1)
                plot2.addLayer(layB2)
                plot2.addLayer(layB3)
                plot2.addLayer(layB4)
                ymin = boxMinMax[0]
        #Loop over pixels
        for j in range (cube.dimensions[1]):
            for k in range (cube.dimensions[2]): 
                specIn = cube.getSpectrum1d(j, k)
                ppTablePos = ppTableIn['HDU_1']['row'].data.where((ppTableIn['HDU_1']['row'].data == j).and(ppTableIn['HDU_1']['column'].data == k))
                if ppTablePos.toInt1d().size < 1:
                    continue
                ppTable = ppTableIn['HDU_1'].select(ppTablePos)
                if array == 'SLW':
                    specIn.meta['channelName'] = StringParameter('SLWC3')
                else:                 
                    specIn.meta['channelName'] = StringParameter('SSWD4')
                #Check for NaNs
                #if specIn.flux.where(IS_FINITE(specIn.flux)).toInt1d().size > 0:
                if specIn.flux.where(IS_NAN(specIn.flux)).toInt1d().size == 0:
                    '''
                    ## Main function is `catalogueFeatures' ##
                    ## Which calls functions in this order: ##
                       # For inital fit and subtraction of the contiuum
                         - findJumps: 
                             Looks for jumps to indicate significant features
                             Flags and weights spectrum for fitting
                         - removeContinuum:
                             Initial fit to the continuum
                             Outputs a continuum subtracted spectrum and the fit
                      # The finding and fitting process begins and loops over the snrThresholds
                         - listFeatures: 
                             finds new features > snrThreshold
                         - fitFeatures: 
                             fits the old+new features and keeps new features > snrThreshold
                      # Outputs are then prepared, including taking the final SNRs
                         - findResidualSnr: 
                             takes the SNR using the fitted peak and local 
                             noise from the residual
                    '''
                    try:
                        features, snrs, snrsIter, cutFeatures, cutSnrs, cutSnrsIter, thresholds, \
                    cutThresholds, finalModelArray, finalResidualArray, fitter, continuumModel,\
                    sincModels, cutSincModels, featuresErr, cutFeaturesErr, initialContinuumArray, \
                    peaks, cutPeaks, cutContinuumValues, outputInitialFeatureFreq, \
                    cutOutputInitialFeatureFreq, indexFeatures, cutIndexFeatures,\
                    threshSincs,threshResiduals,threshTotalModels,\
                    gOfFthreshold, linesForGofFthresholds, allSnrForGofFthresholds,\
                    flags, cutSnrs2, cutSnrs3,\
                         residualDet,residualDet2,residualDet3 = \
                              catalogueFeaturesWithPp(specIn,snrThresholds,signs,jumpThreshold,array,flags,\
                                          resampleResolution=resampleResolution,flagWidths=flagWidths, \
                                          polyOrder=polyOrder, snrCut=snrCut, snrRange=snrRange,\
                                          mergeRadius=mergeRadius, subtractBaseline=True, \
                                          baselineOrder=baselineOrder, checkCutResults = checkCutResults,\
                                          avoid=avoid, snrCap=snrCap, testGofFcut=testGofFcut,\
                                          testCheckFit=testCheckFit,testNeutralC=testNeutralC,\
                                          estimateRedshift=estimateRedshift,testGofFcutLines=testGofFcutLines,\
                                          maxIter = maxIter,limitValue=limitValue,checkValue=checkValue,\
                                          apod=apod, noExtraBox=False, ppTable = ppTable)
                    except:
                        continue
                    if checkCutResults == False:
                        #Append the detector array to the feature index
                        for ll in range(indexFeatures.size):
                            indexFeatures[ll]=indexFeatures[ll]+array+'_%s_%s'%(j,k)
                        for ll in range(cutIndexFeatures.size):
                            cutIndexFeatures[ll]=cutIndexFeatures[ll]+array+'_%s_%s'%(j,k)     
                    #Get the continuum from the results
                    param = continuumModel.getFittedParameters()
                    paramError = continuumModel.getFittedStdDev()
                    poly = PolynomialModel(polyOrder)
                    poly.setParameters(param)
                    con1 = poly(specIn.wave)
                    continuum.image[:,j,k] = con1
                    continuumCat['array'].append(array)
                    continuumCat['row'].append(j)
                    continuumCat['column'].append(k)
                    continuumCat['ra'].append(specIn.meta['ra'].value)
                    continuumCat['dec'].append(specIn.meta['dec'].value)
                    continuumCat['param'].append(param)
                    continuumCat['error'].append(paramError)
                    #for addPara in range(param.size):
                    #    continuum.meta['%s_%s_p%s'%(j,k,addPara)]=DoubleParameter(param[addPara],'Coefficient p%s of order %s polynomial fitted to the [%s,%s] spaxel continuum'%(addPara,polyOrder,j,k))
                    ###GofF per pixel
                    if testGofFfinal:
                        ###When using the final fitted features, GoF can be run outside the main function
                        specInSub = specIn.copy()
                        specInSub.flux = specInSub.flux-con1
                        gOfFpix = gOfF(cutFeatures,specInSub,cutSincModels,rng=gOfFRng)
                    ###Line evidence per pixel
                    if testFitEvidence:
                        if limitFitEvidence:
                            if cutFeatures.size > 35:
                                print 'Number of features found for %s %s (%s,%s) is > 35, so skipping ToFE'%(obsid,array,j,k)
                                oddsPix = Double1d(cutFeatures.size)-99.0
                                resultPix = oddsPix.copy()
                            else:
                                resultPix, oddsPix = fitEvidence(specIn,cutFeatures,normalise=True, plotIt=plotFitEvidence, verbose=verboseFitEvidence, useWeights=useWeights, polyOrder=polyOrder)
                        else:
                            resultPix, oddsPix = fitEvidence(specIn,cutFeatures,normalise=True, plotIt=plotFitEvidence, verbose=verboseFitEvidence, useWeights=useWeights, polyOrder=polyOrder)
                    if checkCutResults:
                        print 'Mid loop products are available for snrThreshold %s:'%(checkCutResults)
                        print '    tempTotalMod, tempResidual, cutTotalModels, cutResiduals'
                        print ' **Only cutSnrs is filled at this point (not snrs)**'
                        tempTotalMod = model #checkCutResults total model
                        tempResidual = residual #checkCutResults residual
                        cutTotalModels = cutFeatures #All total models up to checkCutResults threshold
                        cutResiduals = cutSnrs #All residuals up to checkCutResults threshold
                        #Also: model, residual, fit, finalContinuum, models
                        stop
                    #Print the number of sources above the SNR cut were found
                    if printCubeVerbose:
                        print 'SNR cut > %s found %s features for %s (%s,%s)'%(snrCut,cutFeatures.size,array,j,k)
                    nFeatures[getSpecGo]+=cutFeatures.size
                    ###SORT THE CUT RESULTS
                    cutFeatures, [cutIndexFeatures,cutIndexFeatures,cutFeaturesErr,cutSnrs3,\
                        cutSnrsIter, cutThresholds,cutPeaks,cutContinuumValues,cutOutputInitialFeatureFreq], \
                             indOut = sortArrays(cutFeatures,[cutIndexFeatures,cutIndexFeatures,\
                                   cutFeaturesErr,cutSnrs,cutSnrsIter,cutThresholds,\
                                        cutPeaks,cutContinuumValues,cutOutputInitialFeatureFreq])
                    ###Append output for prototype product
                    freqArrayPix+=cutFeatures
                    indArrayPix+=cutIndexFeatures
                    freqErrArrayPix+=cutFeaturesErr
                    initialFreqArrayPix+=cutOutputInitialFeatureFreq
                    snrArrayPix+=cutSnrs
                    snrIterArrayPix+=cutSnrsIter
                    arrayArrayPix+=[array]*cutFeatures.size
                    rowArrayPix+=[j]*cutFeatures.size
                    columnArrayPix+=[k]*cutFeatures.size
                    raArrayPix+=[specIn.meta['ra'].value]*cutFeatures.size
                    decArrayPix+=[specIn.meta['dec'].value]*cutFeatures.size
                    threshArrayPix+=cutThresholds
                    modelsArrayPix+=cutSincModels
                    residArrayPix.append(threshResiduals)
                    totalModArrayPix.append(threshTotalModels)
                    sincArrayPix.append(threshSincs)
                    rowC = rowCs[getSpecGo]
                    colC = colCs[getSpecGo]
                    if (j == rowC) & (k == colC):
                        #zEstimateDets.append(zEstimate*c)
                        featuresCentDets+=cutFeatures
                        errCentDets+=cutFeaturesErr
                        snrCentDets+=cutSnrs
                        #if zEstimate > -99:
                        #    zEstimateAve+=0.5*zEstimate*c
                    ###For the combined output for all obsids
                    freqDetsAll+=cutFeatures
                    indDetsAll+=cutIndexFeatures
                    freqErrDetsAll+=cutFeaturesErr
                    snrDetsAll+=cutSnrs
                    detDetsAll+=[array]*cutFeatures.size
                    raAll+=[specIn.meta['ra'].value]*cutFeatures.size
                    decAll+=[specIn.meta['dec'].value]*cutFeatures.size
                    rowAll+=[j]*cutFeatures.size
                    columnAll+=[k]*cutFeatures.size
                    ###And just for the combined flags
                    ###the features per single array
                    featuresArray+=cutFeatures
                    if testGofFfinal:
                        gOfFpix = gOfFpix[Selection(indOut)]
                        gOfFArrayAll+=gOfFpix
                        gOfFAll+=gOfFpix
                    if testFitEvidence:
                        oddsPix = oddsPix[Selection(indOut)]
                        resultPix = resultPix[Selection(indOut)]
                        oddsArrayAll+=oddsPix
                        resultArrayAll+=resultPix
                        oddsAll+=oddsPix
                        resultAll+=resultPix
                    if plotPostCards:
                        #rowC = cube.wcs.naxis1/2
                        #colC = cube.wcs.naxis2/2 
                        rowC = rowCs[getSpecGo]
                        colC = colCs[getSpecGo]
                        if (j == rowC) & (k == colC):
                            lineThicks = Double1d(cutFeatures.size)+0.75
                            if array == 'SLW':
                                cSpax += '[%s,%s](%s); '%(j,k,cutFeatures.size)
                                linCol = red
                                specCol = darkRed
                                posOverLap = cutFeatures.where(cutFeatures > overlap[0])
                                lineThicks[posOverLap] = 1.5
                            else:
                                cSpax += '[%s,%s](%s)'%(j,k,cutFeatures.size)
                                linCol = tickBlue
                                specCol = darkBlue
                            plot2.addLayer(LayerXY(specIn.wave,specIn.flux,color=specCol))
                            ##add the lines
                            for go in range(cutFeatures.size):
                                lin = cutFeatures[go]
                                snr = cutSnrs[go]
                                lineThick = lineThicks[go]
                                totalPeak = cutPeaks[go] + cutContinuumValues[go]
                                subPeak = cutPeaks[go]
                                if snr > 0:
                                    plot2.addLayer(LayerXY(Double1d([lin,lin]),Double1d([totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05,totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05*LOG(snr)]),\
                                                color=linCol,stroke = lineThick))
                                else:
                                    ymin = boxMinMax[0]-(boxMinMax[1]-boxMinMax[0])*0.05
                                    plot2.addLayer(LayerXY(Double1d([lin,lin]),Double1d([totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05*LOG(ABS(snr)),totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05]),\
                                                color=linCol,stroke = lineThick))
                            plotCon = continuum.image[:,rowC,colC]
                            if array == 'SSW':
                                maxCon = MAX(plotCon)
                            plot2.addLayer(LayerXY(specIn.wave,con1,color=Color.GREEN,stroke=1))
                else:
                    if printCubeVerbose:
                        print 'Skipping (%s,%s)'%(j,k)
        print 'SNR cut > %s found %s features for %s'%(snrCut,nFeatures[getSpecGo],array)
        message+='%s %s features found, '%(nFeatures[getSpecGo],array)
        if saveContinuumFits:
            simpleFitsWriter(continuum,'%s%s_%s_continuum.fits.gz'%(continuumPath,obsid,array),compression='GZIP')
    #GoF flags per cube
    if testGofFfinal:
        gOfFflags = getGofFflags(gOfFArrayAll)
        gOfFflagsArrayAll=gOfFflags
        gOfFflagsAll+=gOfFflags
    #ToFE flags per cube
    if testFitEvidence:
        fitEvidenceFlags = getFitEvidenceFlags(resultArrayAll)
        fitEvidenceFlagsArrayAll=fitEvidenceFlags
        fitEvidenceFlagsAll+=fitEvidenceFlags
    #Combined flags
    if testGofFfinal:
        combinedFlags = getCombinedFlags(gOfFArrayAll, resultArrayAll)
        combinedFlags = getCombinedFlags2(combinedFlags,featuresArray,overlap,noisyRegions)
        combinedFlagsArrayAll=combinedFlags
        combinedFlagsAll += combinedFlags
    #######PREPARE AND SAVE THE OUTPUT PRODUCT PER OBSID#######
    featureCatPerObsid = prepareOutputTable([Double1d(freqArrayPix),Double1d(freqErrArrayPix),Double1d(snrArrayPix),\
              String1d(arrayArrayPix),Long1d(rowArrayPix),Long1d(columnArrayPix),Double1d(raArrayPix),Double1d(decArrayPix)],\
              ['frequency','frequencyError','SNR','array','row','column','ra','dec'],\
              units=[Frequency.GIGAHERTZ,Frequency.GIGAHERTZ,None,None,None,None,None,None],\
              description=['Measured frequency','Error on frequency','Peak/local error',\
              'Bolometer array','Cube row','Cube column','RA of spaxel','Dec of spaxel'],\
              addMeta = True, metas=[cube.meta['obsid'].value,cube.meta['object'].value,\
                           cube.meta['odNumber'].value,cube.meta['processResolution'].value,cube.meta['biasMode'].value,cube.meta['mapSampling'].value,snrCut,avoid,maxCon,nFeatures[0],nFeatures[1],0,0],\
              metaNames=['obsid','object','operationalDay','resolution','biasMode','mapSampling','minimumSnr','edgeMaskWidth','maxContinuum','nFeaturesSlw','nFeaturesSsw'],\
              metaDescription=[cube.meta['obsid'].description,cube.meta['object'].description,cube.meta['odNumber'].description,cube.meta['processResolution'].description,cube.meta['biasMode'].description,\
                           cube.meta['mapSampling'].description,'Minimum SNR limit for features found','Width of spectral edges masked','Maximum of fitted continuum','n features found in SLW','n features found in SSW',\
                           ],\
              metaType=['L','S','L','S','S','S','D','D','D','L','L'],metaUnit=[None,None,None,None,None,None,None,Frequency.GIGAHERTZ,cube.fluxUnit,None,None])
    if testGofFfinal:
        featureCatPerObsid['featureFlag'] = Column(Double1d(combinedFlagsArrayAll))
        featureCatPerObsid['featureFlag'].description = "Feature flag" 
    ###Write out the continuum table if requested
    if saveContinuumTable:
        continuumTable['obsid'] = Column(String1d([obsid]*continuumCat['row'].size))
        continuumTable['array'] = Column(continuumCat['array'])
        continuumTable['row'] = Column(continuumCat['row'])
        continuumTable['column'] = Column(continuumCat['column'])
        continuumTable['ra'] = Column(continuumCat['ra'])
        continuumTable['dec'] = Column(continuumCat['dec'])
        for para in range(polyOrder+1):
            colP = Double1d()
            colE = Double1d()
            for pp in range(len(continuumCat['param'])):
                colP.append(continuumCat['param'][pp][para])
                colE.append(continuumCat['error'][pp][para])
            continuumTable['p%s'%(para)] = Column(colP)
            continuumTable['p%serr'%(para)] = Column(colE)
        asciiTableWriter(file='%s/%s_fittedContinuumParameters.csv'%(continuumPath,obsid),table=continuumTable, writeMetadata=True)      
        simpleFitsWriter(continuum,'%sFits/%s_continuum.fits.gz'%(continuumPath,obsid),compression='GZIP')
    #Estimate redshiftFF
    if len(featuresCentDets) > 0:
        velocity,velocityErr,zFlag,N = redshiftFF(Double1d(featuresCentDets),Double1d(errCentDets),Double1d(snrCentDets))
    else:
        velocity, velocityErr, zFlag, N = -99, -99, 'Nan',-99
    #Adding redshiftFF to individual catalogue metadata
    featureCatPerObsid.meta['radVelEstimate']=DoubleParameter(velocity)
    featureCatPerObsid.meta['radVelEstimate'].description = 'Radial velocity estimate based on features found'
    featureCatPerObsid.meta['radVelEstimate'].unit = Speed.KILOMETERS_PER_SECOND
    featureCatPerObsid.meta['radVelErr']=DoubleParameter(velocityErr)
    featureCatPerObsid.meta['radVelErr'].description = 'Radial velocity error'
    featureCatPerObsid.meta['radVelErr'].unit = Speed.KILOMETERS_PER_SECOND
    featureCatPerObsid.meta['radVelFlag']=StringParameter(zFlag)
    featureCatPerObsid.meta['radVelFlag'].description = 'Radial velocity flag: max number of matches'
    featureCatPerObsid.meta['nMaxMatch']=LongParameter(N)
    featureCatPerObsid.meta['nMaxMatch'].description = 'Number of features with max number of matches'
    #APPEND TO THE FULL OUTPUT
    obsidAll+=[cube.meta['obsid'].value]*len(freqArrayPix)
    resolutionAll+=[cube.meta['processResolution'].value]*len(freqArrayPix)
    biasModeAll+=[cube.meta['biasMode'].value]*len(freqArrayPix)
    mapSamplingAll+=[cube.meta['mapSampling'].value]*len(freqArrayPix)
    #zEstimateAllAve+=[zEstimateAve]*len(freqArrayPix)
    #zEstimateAllSlw+=[zEstimateDets[0]]*len(freqArrayPix)
    #zEstimateAllSsw+=[zEstimateDets[1]]*len(freqArrayPix)
    finalZall+=[velocity]*len(freqArrayPix)
    finalZerrAll+=[velocityErr]*len(freqArrayPix)
    finalZflagAll+=[zFlag]*len(freqArrayPix)
    nMaxMatchAll+=[N]*len(freqArrayPix)
    if writeResultsObsid:
        asciiTableWriter(file='%scatalogues/%s_featuresFound.csv'%(outDir,obsid),table=featureCatPerObsid, writeMetadata=True)
        simpleFitsWriter(featureCatPerObsid,'%scataloguesFits/%s_featuresFound.fits.gz'%(outDir,obsid),compression='GZIP')
    if plotPostCards:
        plot2.setChartFitMode(ChartFitMode.FIT_CONTAINER_SIZE)
        plot2 = addAxis(plot2,xtitle='Frequency [GHz]',ytitle='%s [%s]'%(cube.fluxDescription,cube.fluxUnit.dialogName),\
                 xrange=[430,1580],yrange=[ymin,boxMinMax[1]*0.65],\
                 width=800,height=450,legend=0,legendFS=12,title='%s Observation ID %s %s'%(cube.meta['object'].value,obsid,cSpax))
        plot2.setPlotSize(herschel.ia.gui.plot.renderer.DoubleDimension2D(6, 3))        
        if savePostCards:
            plot2.saveAsPDF('%sPDF/%s_postcard.pdf'%(postCardPlotPath,obsid))
            plot2.saveAsPNG('%s/%s_postcard.png'%(postCardPlotPath,obsid))
        if closePostCards:
            plot2.close()
    if writeLogToFile:
        addToLog(message,logFileName)

if writeLogToFile:
    date = '%s'%(java.util.Date())
    addToLog('-----------------\n%s: Finished\n-----------------'%(date[11:19]),logFileName)

#Create the combined table and write to csv
if len(freqDetsAll):
    featureCatCombined = prepareOutputTable([Double1d(freqDetsAll),Double1d(freqErrDetsAll),Double1d(snrDetsAll),Long1d(obsidAll),String1d(detDetsAll),\
                   Long1d(rowAll),Long1d(columnAll),Double1d(raAll),Double1d(decAll),String1d(biasModeAll),\
                       String1d(resolutionAll),String1d(mapSamplingAll),\
                   Double1d(finalZall),Double1d(finalZerrAll),String1d(finalZflagAll),Int1d(nMaxMatchAll)],\
                   ['frequency','frequencyError','SNR','obsid','array','row','column','ra','dec','biasMode','resolution','mapSampling','radVelEstimate','radVelErr','radVelFlag','nMaxMatch'],\
                   description=['Measured frequency','Error on measured frequency','Peak/local error',cube.meta['obsid'].description,'Bolometer array',\
                         'Spaxel row','Spaxel column','RA','Dec',cube.meta['biasMode'].description,\
                               cube.meta['processResolution'].description,cube.meta['mapSampling'].description,'Radial velocity estimate','Error on radial velocity estimate','Radial velocity flag','n features with max matches'],\
                   units=[Frequency.GIGAHERTZ,Frequency.GIGAHERTZ,None,None,None,None,None,None,None,None,None,None,Speed.KILOMETERS_PER_SECOND,Speed.KILOMETERS_PER_SECOND,None,None],addMeta=True,\
                   metas=[snrCut,avoid],metaNames=['minimumSnr','edgeMaskWidth'],metaUnit=[None,Frequency.GIGAHERTZ],\
                   metaDescription=['Minimum SNR limit for features found','Width of spectral edges masked'],metaType=['D','D'])
    if testGofFfinal:
        featureCatCombined['featureFlag'] = Column(Double1d(combinedFlagsAll))
        featureCatCombined['featureFlag'].description = "Feature flag" 
    if writeResultsCombined:
        nFeatures = featureCatCombined['frequency'].data.size
        asciiTableWriter(file='%scombinedCat_featuresFound_%i_%i_%sobsids.csv'%(outDir,idx0,idxn,len(obsids)),table=featureCatCombined, writeMetadata=True)
else:
    print '\n*No combined catalogue was produced as there were zero features fitted*'
################################################################################
