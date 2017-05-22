################################################################################
####   FTS Spectral Feature Finder functions v24Final 2nd release 30/01/17  ####
################################################################################

######## Functions #######

# 7. Plotting
# 8. Products I/O
# 9.  GoF/LineEvidence
# 10. Check fit/redshift
# 11. Feature finding

################################################################################
####                          6. Import statements                          ####
################################################################################
from java.lang.Math import PI
from herschel.share.unit import *
from java.awt import Color
from herschel.ia.gui.plot.renderer import FillClosureType
from herschel.ia.gui.plot.renderer import PComponentEngine
from herschel.ia.gui.plot.renderer import ChartFitMode
import re
c = Constant.SPEED_OF_LIGHT.value/1000.0

################################################################################
####                              7. Plotting                               ####
################################################################################

###Plotting colours###

teal = Color(0.2,0.4,1.0)
lightTeal = Color(179,226,205)
orange = Color.ORANGE
lightOrange = Color(253,205,172)
red = Color.RED
blue = Color.BLUE
darkRed = Color(0.55,0.0,0.0)
darkBlue = Color(0.0,0.0,0.55)
cyan = Color.CYAN
green = Color.GREEN
darkGreen = Color(0.0,0.6,0.0)
purple = Color(0.8,0.0,1.0)
black = Color.BLACK
lightGray = Color.lightGray
tickBlue = Color(0.0,0.7,1.0)
###Postcard fill colours for the avoid regions
fillColourSlw = Color(0.5,0.0,0.0,0.2)
fillColourSsw = Color(0.0,0.0,0.5,0.2)
###Colours for plotting features found, which are linked to the snrThresholds
threshColours = [black,lightOrange,darkGreen,purple,cyan,teal]
threshColours+=threshColours
#snrThresholds = [100.0, 50.0, 30.0, 10.0, 5.0, 3.0]

#############################################################################
####                           Plotting functions                        ####
#############################################################################

###MAIN DIAGNOSTIC PLOTS
def diagnosticPlots(specIn, conIn, modelIn, features, indexFeatures, snrs, thresholds, plotCol, \
                    sincModels, obsid, name, det_i,\
                    zoomPlots=False,saveMainPlots=False,closePlots=True,plotPath='',\
                    zoomPlotPath=None,checkCutResults=False, snrCut=5.0,\
                    expectedFreq=-99, seePlots=True, res='HR'):
    '''
    Plots the continuum subtracted spectrum, the total model and the residual
    and marks the features foundfeatures found, coloured by < snrCut and snrThreshold
    Input:
        - specIn: spectrum to plot
        - conIn: continuum to subtract (array)
        - modelIn: total fitted model (array)
        - features: frequencies of features found
        - indexFeatures: feature names
        - snrs: SNRs corresponding to features or if checkCutResults these should be snrsIter
        - thresholds: snrThresholds where features were found
        - plotCol: colours linked to the snrThresholds
        - sincModels: fitted sincs
        - obsid: observation ID
        - name: name to use for plot title
        - det_i: detector index, to check if it is the first detector in the list
        - zoomPlots: set to True for zoomPlots
        - saveMainPlots: set to True to save the main plots
        - closePlots: closes the main plots if saveMainPlots==True. Set to False to stop closure
        - plotPath: path where the main plots will be saved if saveMainPlots == True
        - zoomPlotPath: path where the zoom plots will be saved (zoomPlots are always saved)
        - checkCutResults: set to a snrThreshold to obtain the plots at that cut 
        - snrCut: minimum SNR for features  
        - expectedFreq: expected frequency for a doppler shifted CO(7-6) or CO(10-9) feature
        - seePlots: set to False for quicker execution (automatic if closePlots is True)
        - res: for the zoom plot x-range
    Output:
        - plot1: the diagnostic plot   
    '''
    #Create the plot 
    plot1 = PlotXY(seePlots)
    #Mark features < snrCut as dashed grey
    #Mark features > snrCut using colours linked to the snrThresholds
    #See the colour reference plot (black, lightOrange, darkGreen, purple, cyan, teal)
    inLegLow = 1
    inLegHigh = 1
    if expectedFreq > -99:
        plot1.addLayer(LayerXY(Double1d([expectedFreq,expectedFreq]),\
                         Double1d([MIN(specIn.flux-conIn)/3.-3,MAX(specIn.flux-conIn)+1]),\
                         name='Redshift',color=Color.LIGHT_GRAY,stroke=8.0,\
                         line=Style.SOLID,inLegend=inLegLow))
    for ii in range(len(features)):
        feature = features[ii]
        if checkCutResults:
            snr = snrsIter[ii][0]
        else:
            snr = snrs[ii]
        threshold = thresholds[ii]
        col = plotCol[threshold]
        if ABS(snr) < snrCut:
            #Features that didn't make the cut
            plot1.addLayer(LayerXY(Double1d([feature,feature]),\
                             Double1d([MIN(specIn.flux-conIn)/3.-3,MAX(specIn.flux-conIn)+1]),\
                             color=Color.LIGHT_GRAY,stroke=1.25,\
                             line=Style.DASHED,name = 'Features < snrCut',inLegend=inLegLow))
            inLegLow=0
        else:
            #Features included in the output catalogue, coloured by iteration
            plot1.addLayer(LayerXY(Double1d([feature,feature]),\
                             Double1d([MIN(specIn.flux-conIn)/3.-3,MAX(specIn.flux-conIn)+1]),\
                             color=col,stroke=1.25,\
                             name = 'Features >= snrCut',inLegend=inLegHigh))
            inLegHigh = 0
    #plot1.addLayer(LayerXY(specIn.wave,finalResidualArray,name='Residual',color=Color.GREEN))
    plot1.addLayer(LayerXY(specIn.wave,specIn.flux-conIn,name='Spectrum',color=Color.BLUE))
    plot1.addLayer(LayerXY(specIn.wave,modelIn-conIn,name='Fitted model',color=Color.RED))
    ######
    # zoom in on each feature
    ######
    if zoomPlots:
        for go in range(len(features)):
            modPara = sincModels[go].getFittedParameters()
            sinc = SincModel()
            sinc.setParameters(modPara)
            sinc = sinc(specIn.wave)
            layTemp = LayerXY(specIn.wave,sinc,name='Fitted sinc',color=Color.BLACK,line=Style.DASHED,stroke=1.0)
            plot1.addLayer(layTemp)
            feature = features[go]
            featureInd = indexFeatures[go]
            if checkCutResults:
                snr = snrsIter[go][0]
            else:
                snr = snrs[go]
            snrIter = snrsIter[go]
            threshold = thresholds[go]
            if res == 'HR':
                rng = [feature-15,feature+15]
            else:
                rng = [feature-45,feature+45]
            xr=[rng[0],rng[1]]
            sel =  specIn.wave.where((specIn.wave > rng[0]) & (specIn.wave < rng[1]))
            yr=[MIN((specIn.flux-conIn)[sel])-0.1,MAX((specIn.flux-conIn)[sel])+0.2]
            plot1 = addAxis(plot1,xtitle='Frequency [GHz]',ytitle='%s [%s]'%(specIn.fluxDescription,specIn.fluxUnit.dialogName),\
                xrange=xr,yrange=yr,\
                width=0,height=0,legend=1,legendFS=12,title='%s %s %s'%(featureInd, name,det))
            plot1.addAnnotation(Annotation(xr[0]+1.5,yr[1]-0.2*yr[1],'Feature (GHz):',fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+1.5,yr[1]-0.3*yr[1],'SNR:',fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+1.5,yr[1]-0.4*yr[1],'SNR iter:',fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+1.5,yr[1]-0.5*yr[1],'Threshold:',fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+7,yr[1]-0.2*yr[1],'%0.2f'%(feature),fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+7,yr[1]-0.3*yr[1],'%0.1f'%(snr),fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+7,yr[1]-0.4*yr[1],'%0.1f'%(snrIter[0]),fontSize=10))
            plot1.addAnnotation(Annotation(xr[0]+7,yr[1]-0.5*yr[1],'%0.1f'%(threshold),fontSize=10))
            plot1.saveAsPDF('%s%i_%s_%0.1f_%s.pdf'%(zoomPlotPath,obsid,specIn.meta['channelName'].value,feature,featureInd))
            plot1.removeLayer(layTemp)
            plot1.clearAnnotations()
    #rescale plot
    plot1 = addAxis(plot1,xtitle='Frequency [GHz]',ytitle='%s [%s]'%(specIn.fluxDescription,specIn.fluxUnit.dialogName),\
            xrange=[MIN(specIn.wave),MAX(specIn.wave)],yrange=[MIN(specIn.flux-conIn),MAX(specIn.flux-conIn)],\
            width=0,height=0,legend=1,legendFS=12,title='%s %s'%(name,det))
    if saveMainPlots:
        plot1.saveAsPDF('%s%i_%s.pdf'%(plotPath,obsid,det))
        if closePlots:
            plot1.close()
    return plot1

###POSTCARDS
def preparePostcardPlots(spec, avoid, seePostCards=False):
    '''
    Prepares the postcard plot per obsid
    Input:
        - spec: input spectrum
        - avoid: avoided region at the ends of the bands to plot shaded
        - seePostCards: set to False to hide the postcard plots
    Output:
        - plot2: plot to add fitted features and spectrum to
        - boxMinMax: used for setting the yrange
    '''
    plot2 = PlotXY(seePostCards)   
    #Create shaded regions
    boxMinMax = [MIN([MIN(spec[0]['SLWC3'].flux),MIN(spec[0]['SSWD4'].flux)]),\
                  MAX([MAX(spec[0]['SLWC3'].flux)*2,MAX(spec[0]['SSWD4'].flux)*2])]
    boxX = Double1d([0,0,avoid,avoid,0])
    boxY = Double1d([boxMinMax[0]-boxMinMax[1],boxMinMax[1],boxMinMax[1],boxMinMax[0]-boxMinMax[1],boxMinMax[0]-boxMinMax[1]])
    #Add avoid regions, coloured depending on the detector
    layB1 = LayerXY(boxX+spec[0]['SLWC3'].wave[0],boxY,stroke=0,color=Color(0,0,0,0))
    layB1.style.fillClosureType=FillClosureType.SELF
    layB1.style.fillPaint=fillColourSlw
    layB1.style.fillEnabled=1
    layB2 = LayerXY(boxX+spec[0]['SLWC3'].wave[-1]-avoid,boxY,stroke=0,color=Color(0,0,0,0))
    layB2.style.fillClosureType=FillClosureType.SELF
    layB2.style.fillPaint=fillColourSlw
    layB2.style.fillEnabled=1
    layB3 = LayerXY(boxX+spec[0]['SSWD4'].wave[0],boxY,stroke=0,color=Color(0,0,0,0))
    layB3.style.fillClosureType=FillClosureType.SELF
    layB3.style.fillPaint=fillColourSsw
    layB3.style.fillEnabled=1
    layB4 = LayerXY(boxX+spec[0]['SSWD4'].wave[-1]-avoid,boxY,stroke=0,color=Color(0,0,0,0))
    layB4.style.fillClosureType=FillClosureType.SELF
    layB4.style.fillPaint=fillColourSsw
    layB4.style.fillEnabled=1
    plot2.addLayer(layB1)
    plot2.addLayer(layB2)
    plot2.addLayer(layB3)
    plot2.addLayer(layB4)
    return plot2, boxMinMax

def postcardPlots(plot2, spec, det, det_i, cutFeatures, cutSnrs, cutPeaks, \
                  cutContinuumValues, overlap, boxMinMax, withTicks = True):
    '''
    Adds the spectrum to and marks the found features > snrCut on the prepared postcard plot
    Input:
        - plot2: plot prepared by preparePostcarPlots
        - spec: the spectrum to plot
        - det: the detector to plot
        - det_i: to check if its the first detector in the last
        - cutFeatures: frequencies of features found > snrCut
        - cutSnrs: SNRs corresponding to cutFeatures, for scaling the feature marks
        - cutPeaks: fitted peaks corresponding to cutFeatures
        - cutContinuumValues:  continuum values corresponding to cutFeatures
        - overlap: overlap region between the bands
        - boxMinMax: output from preparePostcardPlots and used to define ymin
    '''
    lineThicks = Double1d(cutFeatures.size)+0.75
    if det[:3] == 'SLW':
        linCol = red
        specCol = darkRed
        posOverLap = cutFeatures.where(cutFeatures > overlap[0])
        lineThicks[posOverLap] = 1.5
    else:
        linCol = tickBlue
        specCol = darkBlue
    plot2.addLayer(LayerXY(specIn.wave,spec[0][det].flux,color=specCol))
    ymin = boxMinMax[0]
    ##add the features
    if withTicks:
        for go in range(cutFeatures.size):
            feature = cutFeatures[go]
            snr = cutSnrs[go]
            lineThick = lineThicks[go]
            totalPeak = cutPeaks[go] + cutContinuumValues[go]
            subPeak = cutPeaks[go]
            if snr > 0:
                plot2.addLayer(LayerXY(Double1d([feature,feature]),Double1d([totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05,totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05*LOG(snr)]),\
                            color=linCol,stroke = lineThick))    
                ymin = boxMinMax[0]       
            else:
                ymin = boxMinMax[0]-(boxMinMax[1]-boxMinMax[0])*0.05
                plot2.addLayer(LayerXY(Double1d([feature,feature]),Double1d([totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05*LOG(ABS(snr)),totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05]),\
                            color=linCol,stroke = lineThick))
    if det_i == 0:
        plot2.setChartFitMode(ChartFitMode.FIT_CONTAINER_SIZE)
        plot2 = addAxis(plot2,xtitle='Frequency [GHz]',ytitle='%s [%s]'%(spec[0][det].fluxDescription,spec[0][det].fluxUnit.dialogName),\
                    xrange=[430,1580],yrange=[ymin,boxMinMax[1]*0.65],width=800,height=450,\
                    legend=0,legendFS=12,title='%s Observation ID %s'%(spec.meta['object'].value,obsid))
        plot2.setPlotSize(herschel.ia.gui.plot.renderer.DoubleDimension2D(6, 3))
    return plot2

def postcardPlotsAddSnr(plot2, det, cutFeatures, cutSnrs, cutPeaks, \
                  cutContinuumValues, overlap, boxMinMax, shift=0.0, colourShift=[0.0,0.0,0.0]):
    '''
    Adds the spectrum to and marks the found features > snrCut on the prepared postcard plot
    Input:
        - plot2: plot prepared by preparePostcarPlots
        - spec: the spectrum to plot
        - det: the detector to plot
        - det_i: to check if its the first detector in the last
        - cutFeatures: frequencies of features found > snrCut
        - cutSnrs: SNRs corresponding to cutFeatures, for scaling the feature marks
        - cutPeaks: fitted peaks corresponding to cutFeatures
        - cutContinuumValues:  continuum values corresponding to cutFeatures
        - overlap: overlap region between the bands
        - boxMinMax: output from preparePostcardPlots and used to define ymin
    '''
    lineThicks = Double1d(cutFeatures.size)+0.75
    if det[:3] == 'SLW':
        linCol = Color(0.8-colourShift[0],0.0,0.0)
        specCol = darkRed
        posOverLap = cutFeatures.where(cutFeatures > overlap[0])
        lineThicks[posOverLap] = 1.5
    else:
        linCol = Color(0.0,0.5-colourShift[1],0.8-colourShift[2])
        specCol = darkBlue
    ##add the features
    for go in range(cutFeatures.size):
        feature = cutFeatures[go] + shift
        snr = cutSnrs[go]
        lineThick = lineThicks[go]
        totalPeak = cutPeaks[go] + cutContinuumValues[go]
        subPeak = cutPeaks[go]
        if snr > 0:
            plot2.addLayer(LayerXY(Double1d([feature,feature]),Double1d([totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05,totalPeak+(boxMinMax[1]-boxMinMax[0])*0.05*LOG(snr)]),\
                        color=linCol,stroke = lineThick))    
            ymin = boxMinMax[0]       
        else:
            ymin = boxMinMax[0]-(boxMinMax[1]-boxMinMax[0])*0.05
            plot2.addLayer(LayerXY(Double1d([feature,feature]),Double1d([totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05*LOG(ABS(snr)),totalPeak-(boxMinMax[1]-boxMinMax[0])*0.05]),\
                        color=linCol,stroke = lineThick))
    return plot2


#Other plotting functions
def addAxis(plot,xtitle='Frequency [GHz]',ytitle='Flux density [Jy]',xrange=0,\
            yrange=0,width=0,height=0,legend=0,legendFS=10,title=''):
    '''
    Handy way to add a title, axis titles, axis ranges, 
    resize the plot window and handle the legend
    '''
    plot.xtitle=xtitle
    plot.ytitle=ytitle
    plot.titleText=title
    plot.legend.visible=legend
    plot.legend.fontSize=legendFS
    if xrange:
        plot.xrange=xrange
    if yrange:
        plot.yrange=yrange
    if width:
        plot.width=width
    if height:
        plot.height=height
    return plot

def plotColours(colours,snrThresholds):
    '''
    Plots a colour reference for the snrThreshold colours
    Takes a list of colours and snrThresholds
    '''
    plotCol = PlotXY()
    ystr = ['']
    for go in range(len(colours)):
        plotCol.addLayer(LayerXY(Double1d([0,1]),Double1d([len(colours)-go,len(colours)-go]),stroke=2.0,color=colours[go],name='> %s SNR'%snrThresholds[go]))
        ystr.append('%i'%(snrThresholds[len(colours)-go-1]))
    plotCol = addAxis(plotCol,xtitle='',ytitle='SNR colour',xrange=[0,1],yrange=[0,go+2],\
                    width=0,height=300,legend=1,legendFS=12,title='Colour reference')
    ystr.append('')
    plotCol.xaxis.tick.label.visible=False
    plotCol.xaxis.getTick().setMinorNumber(0)
    plotCol.xaxis.getTick().setNumber(0)
    plotCol.xaxis.getTick().visible=False
    plotCol.yaxis.getTick().setMinorNumber(0)
    plotCol.yaxis.tick.setNumber(go+2)
    plotCol.yaxis.tick.label.setFixedStrings(ystr)
    plotCol.yaxis.tick.label.setOrientation(1)
    return plotCol

################################################################################
####                            8. Products I/O                             ####
################################################################################

###Get product from the archive
def getProduct(obsid, res='HR', useHsa=True, mapping=False, cubeName='', \
                   array=False, useExt=False):
    '''
    Retrives product from the HSA
    Input:
        - obsid: observation ID
        - res: resolution of product
        - useHsa: set to False if the pool is stored locally
        - mapping: set to True to obtain the SLW and SSW cubes, otherwise
                   the pointsource calibrated SDS is extracted
    Output:
        - obs: the observation context
        - product: the extracted product
        - name: name used for plotting
    '''
    product = True
    if useHsa:
        query = MetaQuery(herschel.ia.obs.ObservationContext, 'p', 'p.obsid==%s'%obsid,1)
        refs = hsa_storage.select(query)
        try:
            obs = hsa_storage.load(refs[-1].urn).product
        except:
            product = False
            name = '***Couldnt access 0x%X %i***'%(obsid,obsid)
            return False, product, name
        #obs = getObservation(obsid,useHsa=True)
    else:
        obs = getObservation(obsid)
    if ('FAILED' in obs.quality.meta["state"].string):         
        product = False
        name = "0x%X %i failed, skipping it"%(obsid, obsid)
    if product:
        try:
            if obs.meta['processedAs'].value != 'H+LR':
                if res != obs.meta['processedAs'].value:
                    product = False
                    name = '***Defined resolution for 0x%X %i does not match the observation.***'%(obsid,obsid)
        except:
            product = False
            name = '***No processedAs: failed? 0x%X %i***'%(obsid,obsid)
    if product:
        if mapping:
            #Extract the cube
            #if cubeName == 'cube_convol':
            product = [obs.refs["level2"].product.refs["%s_%s_%s"%(res,'SLW',cubeName)].product]
            product.append(obs.refs["level2"].product.refs["%s_%s_%s"%(res,'SSW',cubeName)].product)
            #else:
                #product = [obs.refs["level2"].product.refs["%s_%s"%(res,'SLW')].product]
                #
                #product.append(obs.refs["level2"].product.refs["%s_%s"%(res,'SSW')].product)
        else:
            if useExt:
                product = obs.refs['level2'].product.refs['%s_spectrum_ext'%(res)].product
            else:
                product = obs.refs['level2'].product.refs['%s_spectrum_point'%(res)].product
        name = '%s %s'%(obs.meta['object'].value,obsid)
    else:
        if mapping:
            product = ['failed']
        else:
            product = 'failed'
    return obs, product, name

def sortArrays(sorter, toBeSorted):
    '''
    Sorts a list of arrays depending on sorter
    Input: 
        - sorter: array to use for sorting
        - toBeSorted: list of arrays to sort following the order of sorter
    Output:
        - sorterSorted: reordered sorter
        - sortedOut: the sorted list of toBeSorted
        - ind: sorting indices
    '''
    ind = SORT.BY_INDEX(sorter.copy())
    sorter = sorter[Selection(ind)]
    sortedOut = []
    for arr in toBeSorted:
        arr = arr[Selection(ind)]
        sortedOut.append(arr)
    return sorter,sortedOut,ind

def prepareOutputTable(columns, names, tableDescription='', description=False, \
                         units=False, addMeta=False, metas=False, metaNames=False,\
                            metaType=False, metaDescription=False, metaUnit=False, addFlag=True):
    '''
    Prepares a TableDataset
    Input:
        - columns: columns to add to the table
        - names: column names
        - description: column description
        - units: column units
        - addMeta: set to True to add metadata using metas, metaNames and metaType
        - metas: list of metadata to add
        - metaNames: list of metadata names, corresponding to metas
        - metaType: list of metadata type, corresponding to metas
                    'S' = StringParameter, 'D' = DoubleParameter, 'L' = LongParameter
        - metaDescription: metas description
        - metaUnit: metas units
     Output:
        - outputTable: filled table
    '''
    outputTable = TableDataset(description = tableDescription)
    if addFlag:
        #Add the flag definitions
        flagName = ['flagGood','flagGoodNoisy','flagPoor','flagPoorNoisy']
        flagDes = ['Good fit in lower noise region','Good fit in noisy region','Poor fit in lower noise region','Poor fit in noisy region']
        flagVal = ['0.0','0.1','1.0','1.1']
        for fName,fDes,fVal in zip(flagName,flagDes,flagVal):
            outputTable.meta['%s'%(fName)]=StringParameter(fVal)
            outputTable.meta['%s'%(fName)].description = fDes
    #
    #if columns:
    for col,nam in zip(columns,names):
        outputTable['%s'%nam] = Column(col)
    if description:
        for des,nam in zip(description,names):
            if des:
                outputTable['%s'%nam].description = des
    if units:
        for unit,nam in zip(units,names):
            if unit:
                outputTable['%s'%nam].unit = unit
    if addMeta:
        for met,metNam,typ in zip(metas,metaNames,metaType):
            if typ == 'S':
                outputTable.meta['%s'%metNam]=StringParameter(met)
            elif typ == 'D':
                outputTable.meta['%s'%metNam]=DoubleParameter(met)
            elif typ == 'L':
                outputTable.meta['%s'%metNam]=LongParameter(met)
        if metaDescription:
            for des,metNam in zip(metaDescription,metaNames):
                if des:
                    outputTable.meta['%s'%metNam].description = des
        if metaUnit:
            for unit,metNam in zip(metaUnit,metaNames):
                if unit:
                    outputTable.meta['%s'%metNam].unit = unit                      
    return outputTable

def printFeaturesVerbose(cutFeatures, cutSnrs, name):
    ''' Prints all features found > snrCut and their finalSnr'''
    print 'In ' + name + ', found:'
    for feature, snr in zip(cutFeatures,cutSnrs):
        print 'feature at '+str(feature)+'GHz, with SNR '+str(snr)
    if len(features) == 0:
        print 'no features'
    return

#for the log file
def addToLog(message,fileName):
    file = open(fileName,'a')
    file.write('%s\n'%message)
    file.close()

################################################################################
####                              9.  GoF/ToFE                              ####
################################################################################

###GoF

def gOfF(featureList, specIn, models, rng=3):
    """
    Calculates the goodness of fit per feature found
    Inputs:
        - featureList: input features found
        - specIn: continuum subtracted spectrum
        - specModel: fitted models
        - rng: +/- range over which the goodness of fit is calculated (GHz) default is 3 GHz
    Returns:
        - gFit: the goodness of fit per feature
    """
    gFit = Double1d()
    #for each feature found, use the individual model
    for i in range(len(featureList)):
        featureFrq=featureList[i]
        specOr = specIn.copy()
        #get somewhere to put the individual feature model
        specModel = specIn.copy()
        #get the fitted sinc for that feature
        modPara = models[i].getFittedParameters()
        sinc = SincModel()
        sinc.setParameters(modPara)
        sinc = sinc(specIn.wave)
        specModel.flux = sinc
        #now make the full model omitting the sinc of the show
        modelSubOne = specIn.copy()
        modelSubOne.flux = modelSubOne.flux-modelSubOne.flux
        for jj in range(len(models)):
            if jj != i:
                modPara = models[jj].getFittedParameters()
                sincTemp = SincModel()
                sincTemp.setParameters(modPara)
                modelSubOne.flux+= sincTemp(specIn.wave)
        #subtract the model omitting the sinc of the show
        specOr.flux = specOr.flux-modelSubOne.flux
        stats1 = statistics(ds=specOr, ranges=(featureFrq-rng, featureFrq+rng))
        stats2 = statistics(ds=specModel, ranges=(featureFrq-rng, featureFrq+rng))
        specTmp = specIn.copy()
        specTmp.flux = ((specOr.flux-stats1.getValueAt(0,0))\
            *(specModel.flux-stats2.getValueAt(0,0))/(stats1.getValueAt(0,1)*stats2.getValueAt(0,1)))
        stats2 = statistics(ds=specTmp, ranges=(featureFrq-rng, featureFrq+rng))
        gFit.append(stats2.getValueAt(0,0))
    return gFit

#Tested. Not found useful
def gOfFbadRemover(gOfFmidCut,featuresFound,featuresFoundErr,continuumPeakResult,\
                     continuumResult,peaksAll,snrIterResult,snrThresholdsResult,\
                     flagWidthsAll,sincModels):
    '''
    Removes features found that are below the goodness threshold
    '''
    pos = gOfFmidCut.where(gOfFmidCut > 0.54)
    if pos.toInt1d().size > 0:
        gOfFmidCut,featuresFound,featuresFoundErr,continuumPeakResult,\
                     continuumResult,peaksAll,snrIterResult,snrThresholdsResult,\
                     flagWidthsAll,sincModels = \
                     gOfFmidCut[pos],featuresFound[pos],featuresFoundErr[pos],continuumPeakResult[pos],\
                     continuumResult[pos],peaksAll[pos],snrIterResult[pos],snrThresholdsResult[pos],\
                     flagWidthsAll[pos],sincModels[pos]
    return gOfFmidCut,featuresFound,featuresFoundErr,continuumPeakResult,\
                     continuumResult,peaksAll,snrIterResult,snrThresholdsResult,\
                     flagWidthsAll,sincModels
 
#def getGofFflags(gOfFdets):
#    '''
#    obtain the GoF flags
#    '''
#    gOfFflagsOut=Int1d()
#    for gg in gOfFdets:
#        if gg >= 0.9:
#            gOfFflagsOut.append(0)
#        elif (gg >= 0.8)&(gg < 0.9):
#            gOfFflagsOut.append(1)
#        elif (gg < 0.8):
#            gOfFflagsOut.append(2)
#    return gOfFflagsOut

def getGofFflags(gOfFdets):
    '''
    obtain the GoF flags
    '''
    gOfFflagsOut=Int1d()
    for gg in gOfFdets:
        if gg >= 0.64:
            gOfFflagsOut.append(0)
        elif (gg < 0.64):
            gOfFflagsOut.append(1)
    return gOfFflagsOut

###Total Fit Evidence

def fitEvidence(spec,featureList,normalise=True, plotIt=True, verbose=True, useWeights=False, polyOrder=3):
    """
     INPUTS:
         - spec: a detector spectrum
         - featureList: list of features found [GHz]
         - normalise: default = True, to normalise the frequency axis to numbers close to 1
         - plotIt: if diagnostic plots are to be produced. The plot shows the data, 
                   the model and the residual. On the residual curve you will see:
                   - filled black triangles for features with high odds of being real
                   - filled red triangles are features that are less evident
                   - open red triangles are dubious
         - verbose: whether to print full details to the console
     OUTPUTS:
         - diff: the difference of the evidence for a model with all features and a 
                 model with the feature in question removed, 
                 i.e. diff[k] is the difference of the evidences for the total model 
                 and the model with k-th feature removed.
         - oddsOut: feature odds
    """
    pi = Math.PI
    oddsOut = Double1d()
    # normalise x to be in a good interval around 1.0
    x = spec.wave
    if (normalise):
        scaleF = 1.0/MAX(x)
    else:
        scaleF = 1.0
    x = x*scaleF # rescale to be numbers close to 1
    y = spec.flux
    yerr = spec['error'].data
    if (useWeights):
        w = 1.0/yerr**2 # weights
    else:
        w = Double1d(yerr.length(),1)
    lineX = featureList.copy()*scaleF
    # Select only spectral features within the data range for this detector
    sel = lineX.where((lineX > MIN(x)) & (lineX < MAX(x)))
    lines = lineX[sel]
    nlines = lines.length()
    # Define the models to fit the data 
    # First, start with a n'th order polynomial for the continuum
    m  = PolynomialModel(polyOrder)
    models = m
    # Get the actual spectral resolution, will need it to fix the sinc line
    spectralResolution = spec.meta['actualResolution'].value*scaleF
    # Create a Sinc model for each spectral feature to be fitted
    for line in lines:
        # Create an initial guess for the model fit parameters
        selLine  = x.where(ABS(x - line) == MIN(ABS(x - line)))
        initalAmplitude  = y[selLine][0]
        # Initial Sinc parameters: [peak, feature frequency, width]
        sincParamGuess = Double1d([initalAmplitude, line, spectralResolution/pi])
        m = SincModel()
        m.setParameters(sincParamGuess)
        # Fix the Sinc width
        m.keepFixed(Int1d([2]))
        models.addModel(m)
    # Carry out the global fit
    fitter = LevenbergMarquardtFitter(x,models)
    param = fitter.fit(y,w)
    totalModel = models(x)
    residual   = y - totalModel
    scale = fitter.autoScale()
    np = models.getNumberOfParameters()
    prior = Double1d( np+1 ) + 1000.0
    logpr = fitter.getEvidence( prior )
    if (plotIt):
        plt = PlotXY()
        l1 = LayerXY(x/scaleF,y)
        l2 = LayerXY(x/scaleF,totalModel)
        l3 = LayerXY(x/scaleF,residual)
        plt.addLayer(l1)
        plt.addLayer(l2)
        plt.addLayer(l3)
        plt.xaxis.title.text = "Frequency (GHz)"
        plt.yaxis.title.text = "Flux (%s)"%spec.getFluxUnit().toString().split('[')[0]
    ### Loop over each feature, removing it and refitting to get the evidence of this new model
    diff = Double1d(nlines)
    if verbose:
        print "Input feature, total model evidence, evidence without the feature, difference, odds"
    for k in range(nlines):
        # remove the k-th feature from the list
        mask = Bool1d(nlines,1)
        mask.set(k,0)
        lineOne = lines.get(mask)
        # add the models 
        m = PolynomialModel(polyOrder)
        models = m
        for line in lineOne:
            # Create an initial guess for the model fit parameters
            selLine  = x.where(ABS(x - line) == MIN(ABS(x - line)))
            initalAmplitude  = y[selLine][0]
            # Initial Sinc parameters: [peak, feature frequency, width]
            sincParamGuess = Double1d([initalAmplitude, line, spectralResolution/pi])
            m = SincModel()
            m.setParameters(sincParamGuess)
            # Fix the Sinc width
            m.keepFixed(Int1d([2]))
            models.addModel(m)
        # Carry out the global fit
        fitter = LevenbergMarquardtFitter(x,models)
        param = fitter.fit(y,w)
        totalModel = models(x)
        residual   = y - totalModel
        scale = fitter.autoScale()
        np0 = models.getNumberOfParameters()
        prior0 = Double1d( np0+1 ) + 1000.0
        logpr0 = fitter.getEvidence( prior0 )
        diff[k] = (logpr - logpr0)
        odds = 10**(logpr - logpr0)
        #print "Freq: %f GHz, total model evidence: %f, lines[k]/scaleF, logpr, np, logpr0, np0, logpr - logpr0, odds
        if (verbose):
            print lines[k]/scaleF, logpr, logpr0, diff[k], odds
        oddsOut.append(odds)
    if (plotIt):
        sel = diff.where(diff >= 1)
        nsel = sel.length()
        if (nsel > 0):
            lix = lines[sel]/scaleF
            l4 = LayerXY(lix,Float1d(nsel))
            l4.setLine(Style.NONE)
            l4.setSymbol(Style.FTRIANGLE)
            l4.setColor(java.awt.Color.BLACK)
            plt.addLayer(l4)
        sel = diff.where((diff < 1).and(diff >= -5))
        nsel = sel.length()
        if (nsel > 0):
            lix = lines[sel]/scaleF
            l5 = LayerXY(lix,Float1d(nsel))
            l5.setLine(Style.NONE)
            l5.setSymbol(Style.FTRIANGLE)
            l5.setColor(java.awt.Color.RED)
            plt.addLayer(l5)
        sel = diff.where((diff < -5))
        nsel = sel.length()
        if (nsel > 0):
            lix = lines[sel]/scaleF
            l6 = LayerXY(lix,Float1d(nsel))
            l6.setLine(Style.NONE)
            l6.setSymbol(Style.TRIANGLE)
            l6.setColor(java.awt.Color.RED)
            plt.addLayer(l6)
    return diff, oddsOut

#def getFitEvidenceFlags(result):
#    '''
#    obtain the total fit evidence flags for diff
#    '''
#    flagsOut = Int1d()
#    for diff in result:
#        if diff >=1:
#            flagsOut.append(0)
#        elif (diff < 1)&(diff >= -5):
#            flagsOut.append(1)
#        elif diff < -5:
#            flagsOut.append(2)
#    return flagsOut

def getFitEvidenceFlags(result):
    '''
    obtain the total fit evidence flags for diff
    '''
    flagsOut = Int1d()
    for diff in result:
        if diff >=-6:
            flagsOut.append(0)
        elif (diff < -6):
            flagsOut.append(1)
    return flagsOut

def getCombinedFlags(gOfFdets, results):
    '''
    Assign a combined GofF and ToFE flag
    '''
    if results == None:
        results = Double1d(len(gOfFdets))
    if len(gOfFdets)!=len(results):
        print 'Length of GofF and ToFE results must be equal'
        stop
    flagsOut = Double1d()
    for go in range(len(gOfFdets)):
        gg = gOfFdets[go]
        diff = results[go]
        if (gg > 0.65) and (diff >= -6):
            flagsOut.append(0)
        else:
            flagsOut.append(1)
    return flagsOut

def getCombinedFlags2(combinedFlags,featuresFound,overlap,noisyRegions):
    '''
    Add to the combined flag depending on the frequency of the feature found
    0.0: good fit in lower noise region
    0.1: good fit in noisy region
    1.0: poor fit in lower noise region
    1.1: poor fit in noisy region
    '''
    #Find positions for noisy regions and add 0.1 to the flags
    addToFlag = Double1d(combinedFlags.size)+0.1
    featuresFound = Double1d(featuresFound)
    pos = featuresFound.where((featuresFound > noisyRegions[0]) & (featuresFound < overlap[0]))
    if pos.toInt1d().size > 0:
        addToFlag[pos]=0.0
    pos = featuresFound.where((featuresFound < noisyRegions[1]) & (featuresFound > overlap[1]))
    if pos.toInt1d().size > 0:
        addToFlag[pos]=0.0        
    return combinedFlags+addToFlag

################################################################################
####                         10. Check fit/redshift                         ####
################################################################################

def determineFreqRange(specIn, featureList, peaks, sample=5.0):
    """
    Determines the frequency range over which to search for a spectral max.
    Inputs:
        - featureList: list of feature centers [GHz]
        - peaks: fitted peaks [flux units]
        - sample: default sample range, which may be modified
    Returns:
        - maxFreq: max frequency of the sample range for a given feature
        - minFreq: min frequency of the sample range for a given feature
    """
    minFreq = []
    maxFreq = []
    peaks = ABS(peaks)
    # min frequencies
    for li in range(len(featureList)):
        if li == 0:
            if featureList[li]-sample < specIn.wave[0]:
                minFreq.append(specIn.wave[0])
            else:
                minFreq.append(featureList[li]-sample)
            continue
        # takes an average between adjacent point weighted by their amplitudes
        if ABS(featureList[li]-featureList[li-1]) <= sample:
            minFreq.append((peaks[li]*featureList[li] + peaks[li-1]*featureList[li-1])/(peaks[li]+peaks[li-1]))
        else:
            minFreq.append(featureList[li]-sample)
    # max frequencies
    for li in range(len(featureList)):
        if li == len(featureList)-1:
            if featureList[li]+sample > specIn.wave[-1]:
                maxFreq.append(specIn.wave[-1])
            else:
                maxFreq.append(featureList[li]+sample)
            continue
        # takes an average between adjacent point weighted by their amplitudes
        if ABS(featureList[li]-featureList[li+1]) <= sample:
            maxFreq.append((peaks[li]*featureList[li] + peaks[li+1]*featureList[li+1])/(peaks[li]+peaks[li+1]))
        else:
            maxFreq.append(featureList[li]+sample)  
    return Double1d(maxFreq),Double1d(minFreq)   

def checkFit(specIn,newModels,fitter,residual,sampleRange,engine=1,testNeutralC=False,\
             expectedFreq=-99,neighbour=False,\
             additMods=[],continuum=None,snrThreshold=3.0, flagWidth = 5.0,\
             noiseFloor=0.0,mapping=False):
    '''
    Checks the fit and attempts to improve it by searching for local amplitude extrema and performing 
    a second fit using this amplitude and the frequency it corresponds to as the initial guesses in the 
    secondary fit.
    Input:
        - specIn: input spectrum
        - newModels: the fitted models for newFeatures
        - fitter: the fitter
        - residual: current residual
        - sampleRange: the range used to search for the spectral maximum
        - engine: fitting engine. Default is LM
        - testNeutralC: set to True to look for CI
        - expectedFreq: expected frequency for a doppler shifted CO(7-6) or CO(10-9) feature
        - neighbour: the frquency difference between it and its troublesome neighbour
        - additMods: fitted models, which may be updated by checkFit
        - continuum: fitted continuum as array
        - snrThreshold: SNR above which features are retained
        - flagWidth: width for masking
        - noiseFloor: used for calculating a dynamic flag width
    Output:
        - residual: residual on subtraction of the total fitted model
        - fitter: the fitter
        - newFeaturesResult: new fitted feature frequencies > snrThreshold [GHz]
        - newSnrResult: new fitted feature SNRs
        - newPeakResult: new fitted peaks [flux units]
        - newModelResult: new fitted models
        - newErrResult: error on new fitted feature frequencies [GHz]
        - newContinuumResult: continuum value at position of new fitted features [flux units]
        - additMods: fitted models, which may have been updated by checkFit
        - flagWidth: width for masking
    '''
    flagWidthT, modelsToAdd, modelsToRemove = [], [], []
    newFeaturesResult ,newErrResult, newSnrResult = [],[],[]
    newPeakResult, newModelResult ,newContinuumResult = [],[],[]
    checkPeaks ,checkFeatures = Double1d(), Double1d()
    polyPar = fitter.getModelsAsProduct()['Model1']['Parameters'].data
    polyT = PolynomialModel(len(polyPar)-2)
    polyT.setParameters(polyPar[:-1])
    polyT = polyT(specIn.wave)
    fixed = [2]
    # get a list of new peaks
    for mm in newModels:
        checkPeaks.append(mm.getFittedParameters()[0])
        checkFeatures.append(mm.getFittedParameters()[1])
    # determine frequency points of the sampling range
    maxCut, minCut = determineFreqRange(specIn,checkFeatures,checkPeaks,sampleRange)
    for cgi in range(len(newModels)):
        m = newModels[cgi]
        params = m.getFittedParameters()
        p1 = params[1]
        p0 = params[0]
        pa = params[1]
        p2 = params[2]
        parametersStdDev = m.getFittedStdDev()
        measuredFreqErr = parametersStdDev[1]
        #create a sinc using the fitted parameters
        sinc = SincModel()
        sinc.setParameters(params)
        sinc = sinc(specIn.wave)
        #find the snr
        snrFull = sinc/specIn['error'].data
        pos = specIn.wave.where(ABS(sinc) == MAX(ABS(sinc)))
        snr = snrFull[pos]
        contVal = continuum[pos]
        if ABS(snr[0]) < maxFitThresh:
            additMods.append(m)
            sincProfile = ABS(p0*p2*1./(specIn.wave-p1))
            try:
                wind = ABS(specIn.wave[sincProfile.where(sincProfile >= noiseFloor).toInt1d()[0]] - p1) #*
            except:
                wind = flagWidth
            flagWidthT.append(wind)
            newFeaturesResult.append(p1)
            newErrResult.append(measuredFreqErr)
            newSnrResult.append(snr)
            newPeakResult.append(p0)
            newModelResult.append(m)
            newContinuumResult.append(contVal)
            continue
        #determine the indecies of the sample range
        maxInd = specIn.wave.where(MIN(ABS(specIn.wave - maxCut[cgi])) == ABS(specIn.wave - maxCut[cgi])).toInt1d()[0]
        minInd = specIn.wave.where(MIN(ABS(specIn.wave - minCut[cgi])) == ABS(specIn.wave - minCut[cgi])).toInt1d()[0]
        specT = specIn.copy()
        # determine the max or min spectral value in sample range
        if testMaxMinAmp:
            maxAmps = [MAX((specIn.flux-polyT)[minInd:maxInd]),MIN((specIn.flux-polyT)[minInd:maxInd])]
        else:
            if p0 > 0.0:
                maxAmps = [MAX((specIn.flux-polyT)[minInd:maxInd])]
            if p0 < 0.0:
                maxAmps = [MIN((specIn.flux-polyT)[minInd:maxInd])]
        residualMaxMin = Double1d([])
        residualOriginal = Double1d([])
        modelsMaxMin = [] 
        for maxAmp in maxAmps:
            posM = specIn.flux.where((specIn.flux-polyT) == maxAmp).toInt1d()[0]
            lineM = specIn.wave[posM]
            rangeT = int(round(window/(specIn.wave[1]-specIn.wave[0])))
            #determine indecies of spectral segment
            minIndex = pos.toInt1d()[0] -rangeT
            maxIndex = pos.toInt1d()[0] +rangeT
            if minIndex < 0.0:
                minIndex = 0
            if maxIndex > (len(specIn.wave)-1):
                maxIndex = len(specIn.wave)-1
            #set up the spectral segment to be fit
            specT = specIn.copy()
            specT = extract(selection=[0], ds=specT, ranges=(specIn.wave[minIndex],specIn.wave[maxIndex-1]))
            #first model residual
            residualOriginal.append(SUM(ABS(residual[minIndex:maxIndex])))
            #new spetrum fitter for sampled spectral max guess
            fitterMax = SpectrumFitter(specT, False)
            fitterMax.useFitter(engine)
            fitterMax.addModel('Polynomial',[3],list(polyPar[:-1]))
            additionalIndex = checkFeatures.where(ABS(checkFeatures-p1) < window)
            additionalCent = checkFeatures[additionalIndex]
            additionalAmp = checkPeaks[additionalIndex] 
            modNew = [None]*len(additionalCent)
            modNew[0] = fitterMax.addModel('sinc', [maxAmp,lineM,params[2]])
            modNew[0].setFixed(fixed)
            #if additional features are present in spectral segment (there are issues with this)
            aii = 1
            for addInd in range(len(additionalCent)):
                if additionalCent[addInd] == pa:
                    continue
                modNew[aii] = fitterMax.addModel('sinc', [additionalAmp[addInd],additionalCent[addInd],params[2]])
                #modNew[aii].setFixed(fixed)
                modNew[aii].setFixed([1,2])
                aii += 1
            fitterMax.doGlobalFit()
            fitterMax.residual()
            #residual for model based off sampled spectral max guess
            residualMaxFit = SUM(ABS(fitterMax.getResidual()))
            residualMaxMin.append(residualMaxFit)
            modelsMaxMin.append(modNew[0])
        indexMaxMin = residualMaxMin.where(MIN(residualMaxMin) == residualMaxMin).toInt1d()[0]
        residualCurrent = residualOriginal[indexMaxMin]
        residualMaxFit = MIN(residualMaxMin)
        bestMaxMinModel = modelsMaxMin[indexMaxMin]
        passMax = False
        if residualMaxFit < residualCurrent:  
                passMax = True         
                residualCurrent = residualMaxFit
                fitter.removeModel(m)
                m = bestMaxMinModel
                params = m.getFittedParameters()
                p0 = params[0]
                p1 = params[1]# get the error on the fitted parameters
                p2 = params[2]
                parametersStdDev = m.getFittedStdDev()
                measuredFreqErr = parametersStdDev[1]
                sinc = SincModel()
                sinc.setParameters(params)
                sinc = sinc(specIn.wave)
                snrFull = sinc/specIn['error'].data
                pos = specIn.wave.where(ABS(sinc) == MAX(ABS(sinc)))
                snr = snrFull[pos]
                contVal = continuum[pos]
        if ABS(snr[0]) >= snrThreshold:
            if testNeutralC and (ABS(p1-expectedFreq) <= 4.0):
                #residualCurrent = residualMaxFit
                residualT = [residualCurrent,residualMaxFit]
                param_plus = [p0*0.3,p1+neighbour,params[2]]
                param_minus = [p0*0.3,p1-neighbour,params[2]]
                paramTs = [list(params),param_plus,param_minus]
                modelsTT = [[m],[modNew[0]]]
                for tt in range(len(paramTs)): 
                    if tt == 0:
                        continue
                    modelsT = [None]*(len(additionalCent)+1)
                    fitterT = SpectrumFitter(specT, False)
                    fitterT.useFitter(engine)
                    fitterT.addModel('Polynomial',[3],list(polyPar[:-1]))
                    modelsT[0] = fitterT.addModel('sinc', paramTs[0])
                    modelsT[0].setFixed(fixed)
                    modelsT[1] = fitterT.addModel('sinc', paramTs[tt])
                    modelsT[1].setFixed(fixed)
                    aii = 2
                    for addInd in range(len(additionalCent)):
                        if additionalCent[addInd] == pa:
                            continue
                        modelsT[aii] = fitterT.addModel('sinc', [additionalAmp[addInd],additionalCent[addInd],p2])
                        modelsT[aii].setFixed(fixed)
                        aii += 1
                    fitterT.doGlobalFit()
                    fitterT.residual()
                    residualT.append(SUM(ABS(fitterT.getResidual())))
                    modelsTT.append(modelsT)
                residualT = Double1d(residualT)
                pick = residualT.where(residualT == MIN(residualT)).toInt1d()[0]
                modelsT = modelsTT[pick]
                if (pick in [1,2]) and (not passMax):
                    fitter.removeModel(m)
                for mT in modelsT:
                    additMods.append(mT)
                    params = mT.getFittedParameters()
                    p0 = params[0]
                    p1 = params[1]# get the error on the fitted parameters
                    p2 = params[2]
                    parametersStdDev = m.getFittedStdDev()
                    measuredFreqErr = parametersStdDev[1]
                    sinc = SincModel()
                    sinc.setParameters(params)
                    sinc = sinc(specIn.wave)
                    snrFull = sinc/specIn['error'].data
                    pos = specIn.wave.where(ABS(sinc) == MAX(ABS(sinc)))
                    snr = snrFull[pos]
                    contVal = continuum[pos]
                    sincProfile = ABS(p0*p2*1./(specIn.wave-p1))
                    try:
                        wind = ABS(specIn.wave[sincProfile.where(sincProfile >= noiseFloor).toInt1d()[0]] - p1)
                    except:
                        wind = flagWidth
                    flagWidthT.append(wind)
                    newFeaturesResult.append(p1)
                    newErrResult.append(measuredFreqErr)
                    newSnrResult.append(snr)
                    newPeakResult.append(p0)
                    newModelResult.append(mT)
                    newContinuumResult.append(contVal)
            else:
                additMods.append(m)
                sincProfile = ABS(p0*p2*1./(specIn.wave-p1))
                try:
                    wind = ABS(specIn.wave[sincProfile.where(sincProfile >= noiseFloor).toInt1d()[0]] - p1)
                except:
                    wind = flagWidth
                flagWidthT.append(wind)
                newFeaturesResult.append(p1)
                newErrResult.append(measuredFreqErr)
                newSnrResult.append(snr)
                newPeakResult.append(p0)
                newModelResult.append(m)
                newContinuumResult.append(contVal)
    NEWModels = [None]*len(additMods)
    fitter = SpectrumFitter(specIn, False)
    fitter.useFitter(engine)
    fitter.addModel('Polynomial',[3],list(polyPar[:-1]))
    iii = -1
    for nM in additMods:
        iii += 1
        params = nM.getFittedParameters()
        p0 = params[0]
        p1 = params[1]
        p2 = params[2]
        parr = [p0, p1, p2]
        NEWModels[iii] = fitter.addModel('sinc', parr)
        NEWModels[iii].setFixed([1,2])
    fitter.doGlobalFit()
    fitter.residual()
    residual = fitter.getResidual()
    return residual, fitter, newFeaturesResult, newSnrResult, newPeakResult,\
           newModelResult, newErrResult, newContinuumResult, additMods, flagWidth

def redshift1(specIn,array):
    """
    Estimates position of CI region by estimating redshift by searching for 
    peaks around the rest frame 12CO frequencies
    Inputs:
        - specIn: The input spectrum (continuum removed)
        - array: Array being assessed
    Returns:
        - expectedFreq: estimated frequency of the average of the 12CO(7-6) and CI transitions
        - neighbour: the frquency difference between 12CO(7-6) and potential CI emission
    """
    if array[:3] == 'SLW':
        step = ABS(specIn.wave[1] - specIn.wave[0])
        maxLines = []
        for i in range(len(SLWFreqs)):
            restF = SLWFreqs[i]
            indF = specIn.wave.where(MIN(ABS(specIn.wave - restF)) == ABS(specIn.wave - restF)).toInt1d()[0]
            indMax = indF + int(20/step)
            if indMax >= len(specIn.wave):
                indMax = len(specIn.wave) - 1
            indMin = indF - int(50/step)
            if indMin < 0:
                indMin = 0
            maxAmp = MAX(specIn.flux[indMin:indMax])
            maxAmpInd = specIn.flux.where(MIN(ABS(specIn.flux - maxAmp)) == ABS(specIn.flux - maxAmp)).toInt1d()[0]
            lineF = specIn.wave[maxAmpInd]
            maxLines.append(lineF)
        #return 808.0 + MEDIAN(806.651806 * (Double1d(maxLines) - SLWFreqs)/Double1d(maxLines)), MEDIAN(Double1d(maxLines)/SLWFreqs)*2.68
        return 808.0* MEDIAN(Double1d(maxLines)/SLWFreqs), MEDIAN(Double1d(maxLines)/SLWFreqs)*2.68, MEDIAN(((SLWFreqs/Double1d(maxLines))-1))
    else:
        return -999, 0.0, -999

def redshift2(specIn,factor=20,mergeRadius=10):
    """
    Estimates position of CI region by estimating redshift. Looks at the relative 
    differences between spectral peaks above a given threshold and compares them to the 
    relative difference of rest 12CO frequencies.
    Inputs:
        - specIn: The input spectrum (continuum removed)
        - factor: factor*spectrum error sets limit above which peaks are searched for
        - mergeRadius: If two peaks are found within this radius, only the most prominent is used
    Returns:
        - expectedFreq: estimated frequency of the average of the 12CO(7-6) and CI transitions
        - neighbour: the frquency difference between 12CO(7-6) and potential CI emission
    """
    try:
        maxMatch = ''
        mergeRange = mergeRadius*2.0
        indexAbove = specIn.flux.where(specIn.flux > factor*MEAN(specIn['error'].data))
        freqAbove = specIn.wave[indexAbove]
        indexAbove = indexAbove.toInt1d()
        freqMin = freqAbove[0]
        lowerIndex = Int1d([indexAbove[0]])
        upperIndex = Int1d([])
        freqDiff, matchSum = [], []
        freqPeak,amplitude = Double1d([]),Double1d([])
        for i in range(len(freqAbove)):
            freq = freqAbove[i]
            if ABS(freq-freqMin) > mergeRange:
                freqMin = freq
                upperIndex.append(indexAbove[i-1])
                lowerIndex.append(indexAbove[i])
        upperIndex.append(indexAbove[-1])
        for i in range(len(lowerIndex)):
            if lowerIndex[i] == upperIndex[i]:
                freqPeak.append(specIn.wave[lowerIndex[i]])
            else:
                freqPeak.append(specIn.wave[specIn.flux.where(specIn.flux == MAX(specIn.flux[lowerIndex[i]:upperIndex[i]]))])
        for line in freqPeak:
            amplitude.append(specIn.flux[specIn.wave.where(specIn.wave == line)])
        #generate empty output arrays
        freqDiff, restFreq = [], Double1d([])
        velocityIter = [6000]
        for MaxVelocity in velocityIter:
            matchSum = Int1d([])
            Tolerance = (MaxVelocity*refArray)/(MaxVelocity+c)
            Tolinary = {refArray[0]:Tolerance[0],refArray[1]:Tolerance[1],refArray[2]:Tolerance[2],refArray[3]:Tolerance[3],\
                        refArray[4]:Tolerance[4],refArray[5]:Tolerance[5],refArray[6]:Tolerance[6],refArray[7]:Tolerance[7],\
                        refArray[8]:Tolerance[8],refArray[9]:Tolerance[9]}
            i = -1
            for freq in freqPeak:
                i += 1
                matchSum.append(0)
                for sign in [1.0,-1.0]:
                    #generate the difference array for a given feature and all the other features
                    inspectLine = (freqPeak-freq)*sign
                    inspectLine = inspectLine[inspectLine.where(inspectLine >= 0.0)]
                    #inspect each difference array to see how well it matches with the reference array
                    for ref in refArray:
                        if len(inspectLine.where(ABS(inspectLine - ref) <= Tolinary[ref]).toInt1d()) > 0:
                            matchSum[i]+=1
            # Subtract 1 due to double counting of 0.0 diff entry
            matchSum -= 1
            maxMatch = MAX(matchSum)
            maxMatchIndex = matchSum.where(matchSum == maxMatch)
        redshiftFreq = freqPeak[maxMatchIndex]
        for freq in redshiftFreq:
            restIndex = ALLFreqs.where(ABS(ALLFreqs - freq) == MIN(ABS(ALLFreqs - freq)))
            restFreq.append(ALLFreqs[restIndex][0])
        kLine,kRest = Double1d([]),Double1d([])
        for i in range(len(restFreq)):
            rest = restFreq[i]
            if rest in kRest:
                continue
            else:
                kRest.append(rest)
            index = restFreq.where(restFreq == rest)
            amp = amplitude[index]
            lines = redshiftFreq[index]
            maxIndex = amp.where(amp == MAX(amp))
            kLine.append(lines[maxIndex])
        if len(redshiftFreq) < 1:
            raise ValueError
        else:
            #return redshiftFreq,restFreq,maxMatch,len(redshiftFreq)
            return 808.0*MEDIAN(kLine/kRest),2.68*MEDIAN(kLine/kRest), MEDIAN(((kRest/kLine)-1))
    except:
        #print 'An error has occured when estimating the redshift'
        return -999, 0.0, -999

def estimateFinalRedshift1(freqPeak,freqSnr,freqErr,window=50.0):
    """
    Matches up features found to NII or 12CO lines. If no match is found for NII, 
    peaks near 12CO are checked. Windows are set up arround rest frequencies and 
    searched for the most prominent feature. 
    *Used when estimateFinalRedshift2 gives poor results*
    Inputs:
        - freqPeak: found features > SNR cut
        - freqSnr: SNR of found features
        - freqErr: feature position uncertainties
        - window: search window
    Returns:
        - maxLines: lines found in the given windows with the highest snr
        - maxSnrs: SNR corresponding to maxLines
        - maxErr: error corresponding to maxLines
        - freqRest: rest frequencies corresponding to maxLines 
    """
    if len(freqPeak.where(ABS(freqPeak-1461.1314) < 30.0).toInt1d()) > 0:
        indexN = freqPeak.where(ABS(freqPeak-1461.1314) == MIN(ABS(freqPeak-1461.1314))).toInt1d()[0]
        maxLines = Double1d([freqPeak[indexN]])
        maxSnrs = Double1d([freqSnr[indexN]])
        maxErr = Double1d([freqErr[indexN]])
        freqRest = Double1d([1461.1314])
    else:
        #generate empty output arrays
        maxLines,maxSnrs,maxErr,freqRest = Double1d([]),Double1d([]),Double1d([]),Double1d([])
        for freq in ALLFreqs:
            #search symmetric window for features
            inRangeIndex = freqPeak.where(ABS(freq-freqPeak) <= window)
            if len(inRangeIndex.toInt1d()) < 1:
                continue
            checkFreq = freqPeak[inRangeIndex]
            checkSnr = freqSnr[inRangeIndex]
            checkErr = freqErr[inRangeIndex]      
            #might want to add a cut off snr value here
            maxSnrIndex = checkSnr.where(checkSnr == MAX(checkSnr))
            #append line with greatest snr in this window
            maxLines.append(checkFreq[maxSnrIndex][0])
            maxSnrs.append(checkSnr[maxSnrIndex][0])
            maxErr.append(checkErr[maxSnrIndex][0])
            freqRest.append(freq)
    if len(maxLines) < 1:
        raise ValueError
    else:
        return maxLines,maxSnrs,maxErr,freqRest

def estimateFinalRedshift2(freqPeak,freqSnr,freqErr):
    """
    Matches up features found to 12CO lines. Procedure is based on comparing the 
    differences between features found and the differences between 12CO rest frequencies.
    If the result is not satisfactory, estimateFinalRedshift1 is used
    Inputs:
        - freqPeak: found features > SNR cut
        - freqSnr: SNR of found features
        - freqErr: feature position uncertainties
    Returns:
        - redshiftFreq: features matched to 12CO
        - snrs: SNRs of redshiftFreq
        - errors: fitted frequency errors associated to redshiftFreq
        - restFreq: 12CO rest frequencies corresponding to redshiftFreq lines
        - maxMatch: maximum number of features matched
        - len(redshiftFreq): number of features with maxMatch matches to 12CO
    """
    #generate empty output arrays
    freqDiff, restFreq = [], Double1d([])
    for MaxVelocity in velocityIter:
        matchSum = Int1d([])
        Tolerance = (MaxVelocity*refArray)/(MaxVelocity+c)
        Tolinary = {refArray[0]:Tolerance[0],refArray[1]:Tolerance[1],refArray[2]:Tolerance[2],refArray[3]:Tolerance[3],\
                    refArray[4]:Tolerance[4],refArray[5]:Tolerance[5],refArray[6]:Tolerance[6],refArray[7]:Tolerance[7],\
                    refArray[8]:Tolerance[8],refArray[9]:Tolerance[9]}
        i = -1
        for freq in freqPeak:
            i += 1
            matchSum.append(0)
            for sign in [1.0,-1.0]:
                #generate the difference array for a given feature and all the other features
                inspectLine = (freqPeak-freq)*sign
                inspectLine = inspectLine[inspectLine.where(inspectLine >= 0.0)]
                #inspect each difference array to see how well it matches with the reference array
                for ref in refArray:
                    if len(inspectLine.where(ABS(inspectLine - ref) <= Tolinary[ref]).toInt1d()) > 0:
                        matchSum[i]+=1
        # Subtract 1 due to double counting of 0.0 diff entry
        matchSum -= 1
        maxMatch = MAX(matchSum)
        maxMatchIndex = matchSum.where(matchSum == maxMatch)
    redshiftFreq = freqPeak[maxMatchIndex]
    snrs = freqSnr[maxMatchIndex]
    errors = freqErr[maxMatchIndex]
    for freq in redshiftFreq:
        restIndex = ALLFreqs.where(ABS(ALLFreqs - freq) == MIN(ABS(ALLFreqs - freq)))
        restFreq.append(ALLFreqs[restIndex][0])
    if len(redshiftFreq) < 1:
        raise ValueError
    else:
        return redshiftFreq,snrs,errors,restFreq,maxMatch,len(redshiftFreq)

#Main redshift implementation function
def redshiftFF(frequency,freqErr,freqSnr):
    """
    Estimates redshift based on final features found with positive snr > snrCut
    estimateFinalRedshift2 is applied first and if the result is not satifactory
    estimateFinalRedshift1 is applied to the feature with the highest snr.
    Inputs:
        - frequency: found features > SNR cut
        - freqErr: feature position uncertainties
        - freqSnr: SNR of found features
    Returns:
        - estimated redshift
        - error on estimated redshift
        - meth: flag indicating the maximum number of matches between the refArray
          and the difference array generated for each line. "M1" is bad "M10" is best.
          "NII" if nitrogen line is used.
        - N: total number of lines with the maximum number of matches. "1" is bad "10" is best
    """
    #Remove negative 
    posSnrIndex = freqSnr.where(freqSnr > 0.0)
    frequency = frequency[posSnrIndex]
    freqErr = freqErr[posSnrIndex]
    freqSnr = freqSnr[posSnrIndex]
    if len(frequency) < 1:
        #print 'No valid emission lines\n'
        #return -999,-999,'nan',-999  
        return 0, 0,'nan',-1
    #Generate empty output variables
    meth,N = 'nan',0
    try:
        #Try first redshift estimate (based on refArray)
        coLines,coSnrs,coErr,coRest,match,N = estimateFinalRedshift2(frequency,freqSnr,freqErr)
        meth='M%i'%match
    except:
        #print 'There has been an Error in method 1'
        #print obsid
        #return -999,-999,'nan',-999  
        return 0, 0,'nan',-1
    if meth == 'M1' and MAX(freqSnr)>=10:
        try:
            #Not enough CO lines to be reliable. Check for NII
            coLines,coSnrs,coErr,coRest = estimateFinalRedshift1(frequency,freqSnr,freqErr,50.0)
            N = len(coLines)
            meth='M1'
            if coRest[0] == 1461.1314:
                #If NII line is found
                meth = 'NII'
        except:
            #print 'There has been an Error in method 1'
            #return -999,-999,'nan',-999  
            return 0, 0,'nan',-1
    #Determine if there are multiple lines attributed to the same rest frequency.
    #If so, propagate lines with highest snr
    N = len(coLines)
    if N > 1:
        kCoLines,kCoSnrs,kCoErr,kCoRest = Double1d([]),Double1d([]),Double1d([]),Double1d([])
        for dub in ALLFreqs:
            match = coRest.where(coRest == dub)
            if len(match.toInt1d()) > 0:
                keepIndex = coSnrs.where(coSnrs == MAX(coSnrs[match]))
                kCoLines.append(coLines[keepIndex])
                kCoSnrs.append(coSnrs[keepIndex])
                kCoErr.append(coErr[keepIndex])
                kCoRest.append(coRest[keepIndex])
        coLines,coSnrs,coErr,coRest = kCoLines,kCoSnrs,kCoErr,kCoRest
    N = len(coLines)
    velocities = (((coRest/coLines)-1)*c)
    if STDDEV(velocities) > 100.0:
        maxIndex = coSnrs.where(coSnrs == MAX(coSnrs))
        freq,err,snr = coLines[maxIndex],coErr[maxIndex],coSnrs[maxIndex]
        if snr > 10.0:
            try:
                #STDEV of velocities too high to be reliable. Check for NII
                coLines,coSnrs,coErr,coRest = estimateFinalRedshift1(freq,snr,err,50.0)
                N = len(coLines)
                meth='M1'
                if coRest[0] == 1461.1314:
                    #If NII line is found
                    meth = 'NII'
            except:
                #print 'There has been an Error in method 1'
                return 0,0,'nan',-1
    if len(coLines) > 1:
        return MEDIAN(((coRest/coLines)-1)*c), STDDEV(((coRest/coLines)-1)*c),meth,N
    else:
        return MEDIAN(((coRest/coLines)-1)*c), MEDIAN((coRest*c*coErr)/(coLines**2)),meth,N
    
    #print coLines
    #return MEDIAN((coRest*c/coLines)-1), MEDIAN((coRest*c*coErr)/(coLines**2)),meth,N

###########
def getLrSensitivityCurve(specIn,det,nRep):
    if det[:3] == 'SLW':
        wave = Double1d([474.688687,524.688687,574.6886871,624.688687,674.688687,724.688687,\
                              774.688687,824.688687,874.688687,924.688687,974.688687,1024.688687])
        noise = Double1d([0.005589289548116073,0.007050785203565035,0.00538479092743205,\
                    0.0027678183497640467,0.0024452789314772737,0.0021398242148737206,\
                    0.0018364384522615057,0.002720349889066923,0.003758280646095778,\
                    0.004029375703674937,0.0064850201607990175,3.764427510732575E-14])
    else:
        wave = Double1d([969.3462427,1019.3462427,1069.3462427,1119.3462427,1169.3462427,\
                      1219.3462427,1269.3462427,1319.3462427,1369.3462427,1419.3462427,1469.3462427,\
                         1519.3462427,1569.3462427])
        noise = Double1d([0.004307472270066939,0.0040769967997953994,0.004416264419098215,\
        0.0031318220560894597,0.004492314740320095,0.004032574376229342,0.0033361176771723538,\
        0.0038103695254187208,0.0035936107295122337,0.0032458640783488463,0.0042881934331586375,\
        0.004473192049105586,0.001780052852451424])
    scaleFactor = SQRT(6.4*2*nRep/3600.)
    interp = LinearInterpolator(wave,noise,1)
    noiseOut = interp(specIn.wave)
    return noiseOut/scaleFactor

################################################################################
####                         11. Feature finding                            ####
################################################################################

###Used to update the frequencies of already found features after doGlobalFit()
def updateFeatureFrequencies(sincModels):
    """
    Inputs
    - sincModels: list of sinc models fitted to the already found features
    Outputs
    - newFeatureFrequencies: updated frequencies for the already found features
    """
    newFeatureFrequencies = []
    for m in sincModels:
        newFeatureFrequencies.append(m.getFittedParameters()[1])
    return newFeatureFrequencies

###Used to update or recreate feature flagging
def flagSpec1d(specIn, features, flags, flagWidth, recreate = False):
    """
    Inputs
    - specIn: input spectrum1d
    - features: input features as a list
    - flags: existing flags Python dictionary
    - flagWidth: the width of flagging required
    - recreate: set to flag the associated detector from scratch
    Outputs
    - flags: updated feature flagging pyDictionary
    """
    if recreate:
        flags[specIn.meta['channelName'].value]['p1']=Double1d(features)
        flags[specIn.meta['channelName'].value]['fL']=Double1d(features)-flagWidth
        flags[specIn.meta['channelName'].value]['fR']=Double1d(features)+flagWidth
    else:
        flags[specIn.meta['channelName'].value]['p1'].append(Double1d(features))
        flags[specIn.meta['channelName'].value]['fL'].append(Double1d(features)-flagWidth)
        flags[specIn.meta['channelName'].value]['fR'].append(Double1d(features)+flagWidth)    
    return flags

###Used for the preliminary fit to the continuum###
def findJumps(specIn, jumpThreshold, resampleResolution=None, weight=True,\
    flag=True, equalWeight=False):
    """
    Find strong features in the spectrum by looking for rapid changes (jumps)
    Inputs:
        - specIn: input spectrum
        - jumpThreshold: size of jump required (in standard deviations)
        - resampleResolution: resolution used for resampling (GHz)
        - weight: set weight column = 0 for strong features if True
        - flag: set flag column = 1 for strong features if True
        - equalWeight: set initial weight to a constant (1.0)
                      set equalWeight=True for apodized spectra      
    Returns:
        - resampSpec: resampled spectrum
        - specOut: original spectrum that has been flagged and weighted
    """
    if equalWeight:
        size = len(specIn.wave)
        specIn.setError(Double1d([1.0] * size))
    #Pick a method for resampling: 'gaussian', 'trapezoidal', or 'euler'
    resampScheme = 'trapezoidal'
    #Resample the spectrum and find the resampled difference spectrum
    if resampleResolution != None:
        freqBin = resampleResolution / 4.0 # half distance between resampled points
        resampSpec = resample(ds=specIn, density=True, resolution=resampleResolution,\
            scheme=resampScheme, unit='GHz')
        diff = resampSpec.copy()
        diff.setFlux(ABS(resampSpec.flux - SHIFT(resampSpec.flux, 1)))
    else:
        resampSpec = None
        freqBin = (specIn.wave[1] - specIn.wave[0])/2.0
        diff = specIn.copy()
        diff.setFlux(ABS(specIn.flux - SHIFT(specIn.flux, 1)))
    #Find the rms for the difference spectrum
    stats = statistics(ds=diff)
    pos = diff.flux.where(IS_NAN)
    diff.flux[pos] = 0
    rms1 = stats['rms_1'].getData()[0]
    featFreq = diff.wave[diff.flux.where(ABS(diff.flux) > jumpThreshold*rms1)]
    specOut = specIn.copy()
    weightOut = specIn.weight
    flagOut = specIn.flag
    #Loop over jumps and prepare flags and weights
    for i in xrange(featFreq.length()):
        near = specOut.wave.where((specOut.wave > featFreq[i] - 3*freqBin) & \
            (specOut.wave < featFreq[i] + freqBin))
        if weight:
            weightOut[near] = 0
        if flag:
            flagOut[near] = 1
    #set the weights and flags in the input spectrum
    specOut.setAutoConversion(False)
    specOut.setFlag(flagOut)
    specOut.setWeight(weightOut)
    return resampSpec, specOut

###Used for the preliminary subtraction of to the continuum###
def removeContinuum(specIn, polyOrder=3):
    """
    Makes a preliminary fit to the continuum prior to feature finding 
    and returns this model and the continuum subtracted spectrum
    Inputs:
        - specIn: input spectrum ->
                    This should be the weighted and flagged output from findJumps
                    But not the resampled output
        - polyOrder: degree of the polynomial to be fitted
    Returns:
        - specOut: residual flux stored in "flux" column
        - fitter: the fitter
        - continuumModel: the fitted polynomial model
        - totalModel: the fitted polynomial as a Double1d()
    """
    #Pick the fitting engine: : 1 = L-M, 2 = Amoeba, 3 = Linear 
    #                           4 = MP, 5 = Conjugated-Grad
    engine = 1
    #Use the SpectrumFitter
    fitter = SpectrumFitter(specIn, False)
    fitter.useFitter(engine)
    #Add a polynominal model
    continuumModel = fitter.addModel('Polynomial', [polyOrder], Double1d(polyOrder+1).toArray())
    #Set weights to avoid strong features (input should be output from findJumps)
    fitter.setWeight(specIn.weight)
    fitter.doGlobalFit()
    fitter.residual()
    totalModel = fitter.getTotalModel()
    continuumSubtracted = fitter.getResidual()
    #Create the continuum subtracted output
    specOut = specIn.copy()
    specOut.setFlux(continuumSubtracted)
    specOut.setDescription('Continuum subtracted spectrum')
    return specOut, fitter, continuumModel, totalModel

#Used to take the final SNR
def findResidualSnr(residualSpec, featureFq, peak, existingSnr, snrRange=[5.0,25.0], \
                    subtractBaseline=True, baselineOrder=3, ignoreFlag=None, \
                    avoid=10.0, noExtraBox=False):
    """
    Calculates SNR using the fitted feature peak and the RMS of the full residual.
    Inputs:
        - residualSpec: residual spectrum
        - featureFq: frequency of the fitted feature
        - peak: fitted peak
        - snrRange: range either side of feature for noise estimate [GHz]
        - subtractBaseline: extra BLS as default with poly order 3
        - baselineOrder: order of the BLS polynominal
        - ignoreFlag: regions to ignore
    Returns:
        - snr: SNR of (fitted peak/residual)
    """
    residual = residualSpec.flux
    wave = residualSpec.wave
    rms = Double1d(len(wave))
    snr = Double1d(len(wave))
    snrSpec = residualSpec.copy()
    #Find the region for BLS and SNR estimate for the inputted feature
    pos = wave.where((ABS(wave-featureFq) < snrRange[1]).and(ABS(wave-featureFq) \
                            > snrRange[0]).and(\
                            residualSpec.wave < (MAX(residualSpec.wave)-avoid)).and(\
                            residualSpec.wave > (MIN(residualSpec.wave)+avoid)).and(\
                            residualSpec.flag != ignoreFlag))
    #If there are not enough positions, increase the snrRange
    if noExtraBox:
        if subtractBaseline:
            try:
                poly = PolynomialModel(baselineOrder)
                poly.setParameters(Double1d(baselineOrder+1))
                fluxFitter = Fitter(wave[pos], poly)
                baseLine = fluxFitter.fit(residual[pos])
                fluxBaseline = poly(wave[pos])
                rms = QRMS((residual[pos]-fluxBaseLine)-MEAN(residual[pos]-fluxBaseline))
            except:
                rms = QRMS(residual[pos]-MEAN(residual[pos]))      
        else:
            rms = QRMS(residual[pos]-MEAN(residual[pos])) 
        if rms != 0:
            snr = peak / rms
        else:
            snr = 0.0        
    else:
        if pos.toInt1d().size < 17:
            pos = wave.where((ABS(wave-featureFq) < (snrRange[1]+20)).and(ABS(wave-featureFq) \
                                > snrRange[0]).and(\
                                residualSpec.wave < (MAX(residualSpec.wave)-avoid)).and(\
                                residualSpec.wave > (MIN(residualSpec.wave)+avoid)).and(\
                                residualSpec.flag != ignoreFlag))
        #If there are now enough positions, proceed, or set the snr to the existing SNR (if > 5)
        if pos.toInt1d().size >= 17:
            if subtractBaseline:
                try:
                    poly = PolynomialModel(baselineOrder)
                    poly.setParameters(Double1d(baselineOrder+1))
                    fluxFitter = Fitter(wave[pos], poly)
                    baseLine = fluxFitter.fit(residual[pos])
                    fluxBaseline = poly(wave[pos])
                    rms = QRMS((residual[pos]-fluxBaseLine)-MEAN(residual[pos]-fluxBaseline))
                except:
                    rms = QRMS(residual[pos]-MEAN(residual[pos]))      
            else:
                rms = QRMS(residual[pos]-MEAN(residual[pos])) 
            if rms != 0:
                snr = peak / rms
            else:
                snr = 0.0
        else:
            snr = 0.0
    return snr

###Main finding functions###

def listFeatures(specSNR, snrThreshold, flags, mergeRadius=10.0, avoidFlag=1, sign='pos'):
    """
    Compiles a list of potential features found in the SNR spectrum.
    Inputs:
        - specSNR: input SNR spectrum (flux/error column)
        - snrThreshold: SNR threshold for potential feature
        - mergeRadius: width within which features are merged [GHz]
        - avoidFlag: value of flag to avoid
        - sign: indicates +ve or -ve features to find (use 'pos' or 'neg')
    Returns:
        - featureFreq: list of potential new feature frequencies
        - featureSNR: list of potential new feature SNRs
    """
    #Use the SNR spectrum to find potential peaks above the inputted snrThreshold
    if sign == 'pos':
        peaksAll = specSNR.flux.where((specSNR.flux > snrThreshold) & (specSNR.flag != avoidFlag))
    elif sign == 'neg':
        peaksAll = specSNR.flux.where((specSNR.flux < (-1.0*snrThreshold)) & (specSNR.flag != avoidFlag))
    else:
        raise ValueError("listFeatures needs 'sign' setting to 'pos' or 'neg'.")
    #peaks = specSNR.flux.where((ABS(specSNR.flux) > snrThreshold) & (specSNR.flag != avoidFlag))
    #peaksAll = specSNR.flux.where((ABS(specSNR.flux) > snrThreshold) & \
    #                                     (specSNR.flag != avoidFlag)).toInt1d()
    #Do this properly     
#peaksAll = specSNR.flux.where((ABS(specSNR.flux) > snrThreshold).and(\
#                            specSNR.wave < (MAX(specSNR.wave)-avoid)).and(\
#                            specSNR.wave > (MIN(specSNR.wave)+avoid))).toInt1d()
    peaks = Int1d()
    #see if any peaks sit within an [fL:fR] range. 
    for ff in peaksAll.toInt1d():
        pp = specSNR.wave[ff]
        #look for +ve pL and -ve pR
        pL = pp - flags[specSNR.meta['channelName'].value]['fL']
        pR = pp - flags[specSNR.meta['channelName'].value]['fR']
        pos =  pL.where((pL > 0) & (pR < 0))
        if pos.length():
            continue
        else:
            peaks.append(ff)
    peaks = Selection(peaks)
    nPeaks = peaks.length()
    freqPeaks = specSNR.wave[peaks]
    snrPeaks = specSNR.flux[peaks]
    #Sort the output by SNR
    if sign == 'pos':
        sortedIndices = REVERSE(SORT.BY_INDEX((snrPeaks)))
    else:
        sortedIndices = REVERSE(SORT.BY_INDEX(ABS(snrPeaks)))
    sortedFreq = freqPeaks[Selection(sortedIndices)]
    sortedSNR = snrPeaks[Selection(sortedIndices)]
    #snrPeaks = snrPeaks[Selection(sortedIndicies)]
    processed = Bool1d(nPeaks)
    featureFreq = []
    featureSNR = []
    #Loop over the peaks found and merge the fainter nearby peaks
    for i in xrange(nPeaks):
        if processed[i]:
            continue
        nearby = sortedFreq.where(ABS(sortedFreq - sortedFreq[i]) < mergeRadius)
        featureFreq.append(sortedFreq[i])
        featureSNR.append(sortedSNR[i])
        processed[nearby] = True
    if specSNR.meta['obsid'].value == 1342197466:
        if snrThreshold == 50:
            if specSNR.meta['array'].value == 'SLWC3':
                for af in [709,890.8,881,979.3]:
                    featureFreq.append(af)
                    featureSNR.append(51)
        if snrThreshold == 100:
            if specSNR.meta['array'].value == 'SSWD4':
                for af in [1063.07]:
                    featureFreq.append(af)
                    featureSNR.append(110)
    if specSNR.meta['obsid'].value == 1342248242:
        if snrThreshold == 50:
            if specSNR.meta['array'].value == 'SLWC3':
                for af in [805.9]:
                    featureFreq.append(af)
                    featureSNR.append(51)
    if specSNR.meta['obsid'].value == 1342216879:
        if snrThreshold == 100:
            if specSNR.meta['array'].value == 'SLWC3':
                for af in [576.5]:
                    featureFreq.append(af)
                    featureSNR.append(110)
    if specSNR.meta['obsid'].value == 1342193670:
        if snrThreshold == 100:
            if specSNR.meta['array'].value == 'SSWD4':
                for af in [1266.9]:
                    featureFreq.append(af)
                    featureSNR.append(101)
    #Now sort the output by frequency
    sortedIndices = SORT.BY_INDEX(featureFreq)
    featureFreq = list(Double1d(featureFreq)[Selection(sortedIndices)])
    featureSNR = list(Double1d(featureSNR)[Selection(sortedIndices)])
    return featureFreq, featureSNR

def fitFeatures(specIn, residualFlag2, residualFlag3, newFeatures, newSnr, continuum, \
    flags, sincModels, featuresFound, snrThresholdsResult, snrIterResult,\
    fwhm=None,\
    snrThreshold=3.0, flagValue=1, flagWidth=5.0, apod=False, fitter=None,\
    snrCap=500, testCheckFit=False, testNeutralC=False,\
    expectedFreq=-99, neighbour=False, additMods=[], noiseFloor=0.0, \
    avoid=None, limitValue=2.0, maxIter=None, sign=None, checkValue=2.0, ppTable=False):
    """
    Adds new features to the total model and fits to all features
    Inputs:
        - specIn: spectrum to be fitted
        - newFeatures: position of new features to add to model [GHz]
        - newSnr: SNR of new features 
        - continuum: fitted continuum (array)
        - flags: pyDictionary containing the found feature flagging
        - sincModles: a list of all sinc models already fitted
        - fwhm: FWHM of sinc profiles [GHz]. If None then actualResolution is used
        - snrThreshold: minimum fitted SNR for features
        - flagValue: value for flagging at the region on fitted features > snrThreshold
        - flagWidth: HALF the width of region flagged for features found [GHz]
        - apod: set to True if specIn is an apodized spectrum
        - fitter: existing fitter to add the new features to (or to create if None)
        - snrCap: maximum SNR allowed when calculating the initial peak guess
        - testCheckFit: set to True to check and update the fitted peak
        - testNeutralC: set to True to look for the CI feature
        - expectedFreq: expected frequency for a doppler shifted CO(7-6) or CO(10-9) feature
        - neighbour: the frquency difference between it and its troublesome neighbour
        - additMods: fitted models, which may be updated by checkFit
        - noiseFloor: used for calculating a dynamic flag width
        - avoid: avoid width at the edge of the bands [GHz]
        - limitValue: to restrict the initial feature position during the global fit        
        - maxIter: maximum fitting iterations allowed
        - sign: indicates +ve or -ve features to keep (use 'pos' or 'neg')
        - checkValue: if a feature has moved greater than this, remove it
    Returns:
        - specOut: residual spectrum with flagging updated for new features
        - fitter: instance of the spectrumFitter used to fit old+new features and polynomial
        - newFeaturesResult: new fitted feature frequencies > snrThreshold [GHz]
        - newSnrResult: new fitted feature SNRs
        - newPeakResult: new fitted peaks [flux units]
        - newModelResult: new fitted models
        - newErrResult: error on new fitted feature frequencies [GHz]
        - continuumResult: continuum value at position of new fitted features [flux units]
        - additMods: fitted models, which may be updated by checkFit
        - flags: pyDictionary containing the found feature flagging
        - alreadyFoundFeatures: updated freqencies of found features 
    """
    alreadyFoundFeatures = updateFeatureFrequencies(sincModels)
    if avoid == None:
        avoid = 10.0
    #Choose the fitting engine: 1 = L-M, 2 = Amoeba, 3 = Linear
    #                           4 = MP, 5 = Conjugated-Grad
    engine = 1
    N = len(newFeatures)
    newModels = [None] * N
    newFeaturesResult ,newErrResult, newSnrResult = [],[],[]
    newPeakResult, newModelResult ,newContinuumResult = [],[],[]
    fL = flags[specIn.meta['channelName'].value]['fL']
    fR = flags[specIn.meta['channelName'].value]['fR']
    if fwhm == None:
        #deltaSigma = specIn.meta['actualResolution'].value
        deltaSigma = 1.18448225
    else:
        deltaSigma = fwhm / 1.20671
    #Select the fitter if on 1st iteration
    if fitter == None:
        fitter = SpectrumFitter(specIn, False)
        fitter.useFitter(engine)
    #set the fitting tolerance
    if specIn['flux'].unit.name == "Jy":
       tolerance = 1.0E-10
    else:
       tolerance = 1.0E-40
    fitter.setTolerance(tolerance)
    #Set the maximum number of iterations
    if maxIter != None:
        fitter.setMaxIter(maxIter)
    #Loop over the inputted features to add to the total model
    if apod:
        p2 = 1.5*1.20671*deltaSigma/(2.0*SQRT(2.0*LOG(2.0)))
    else:
        p2  = deltaSigma / PI
    for i in xrange(N):
        #DO NOT FIX THIS. SCRIPT FAILS TO FIT FOR LOWER THRESHOLDS IF FIXED
        #THE SCRIPT IS ALSO (INADVERTANTLY) TWEAKED TO TAKE THIS INTO CONSIDERATION LATER ON
        if ppTable:
            p0 = newSnr[i]
        else:
            if newSnr[i] <= 500:
                p0 = newSnr[i] * specIn['error'].data[i]
            else:
                p0 = snrCap * specIn['error'].data[i]
        p1 = newFeatures[i]
        params = [p0, p1, p2]
        #fixed = [2]
        if apod:
            newModels[i] = fitter.addModel('gauss', params)
        else:
            newModels[i] = fitter.addModel('sinc', params)
        newModels[i].setFixed([2])
        #limit the feature positions, to stop them wander off and even falling out
        #of the frequency bands
        if extraConstraints:
            if limitValue != None:
                newModels[i].setLimits(1,p1-limitValue,p1+limitValue)
#        #limit the feature potions to stop them wandering into flagged regions, 
#        #if no nearby flagged regions in one direction, make sure to avoid the 
#        #avoid region
#        #find the nearest fR on the left, or use avoid regions
#        ###NEW BIT limit only those features near to a flagged zone, withing the flagWidth?
#        limitsApplied = False
#        limitOnLimit = 8.0
#        if limitFeature != None:
#            if limitFeature:
#                limitsApplied = True
#                leftLimit = MIN(specIn.wave)+avoid
#                if fR.size > 0:
#                    try:
#                        posLeftLimit = fR.where(min(iiii for iiii in (p1-fR) if iiii > 0) == (p1-fR))
#                        if p1-fR[posLeftLimit][0] > limitOnLimit:
#                            leftLimit = fR[posLeftLimit][0]
#                            #print 'left limit applied'
#                    except:
#                        pass
#                #find nearest fL on the right, or use avoid
#                rightLimit = MAX(specIn.wave)-avoid
#                if fL.size > 0:
#                    try:
#                        posRightLimit = fL.where(min(iiii for iiii in (fL-p1) if iiii > 0) == (fL-p1))
#                        if fL[posRightLimit][0]-p1 > limitOnLimit:
#                            rightLimit = fL[posRightLimit][0]
#                            #print 'right limit applied'
#                    except:
#                        pass       
#                newModels[i].setLimits(1,leftLimit,rightLimit)
#    #Perform the fit
#    removeLimits = 0
    try:
        fitter.doGlobalFit()
    except:
#        #Try removing the limits, if set
#        if limitsApplied:
#            for mgo in range(len(newModels)):
#                newModels[mgo].setLimits(1,0,0)
#            try:
#                fitter.doGlobalFit()
#                removeLimits = 1
#            except:
#                for mgo in range(len(newModels)):
#                    m = newModels[mgo]
#                    fitter.removeModel(m)
#                    featuresOut = updateFeatureFrequencies(sincModels)
#                return specIn, fitter, newFeaturesResult, newSnrResult, newPeakResult, newModelResult, \
#                   newErrResult, newContinuumResult, additMods,\
#                   flags, featuresOut, residualFlag2, residualFlag3,\
#                   removeLimits
#        else:
        for mgo in range(len(newModels)):
            m = newModels[mgo]
            fitter.removeModel(m)
            featuresOut = updateFeatureFrequencies(sincModels)
        return specIn, fitter, newFeaturesResult, newSnrResult, newPeakResult, newModelResult, \
           newErrResult, newContinuumResult, additMods,\
           flags, featuresOut, residualFlag2, residualFlag3
    fitter.residual()
    #Find the total model and residual
    totalModel = fitter.getTotalModel()
    residual = fitter.getResidual()
    #Loop over the fitted newModels to get the fitted parameters 
    #the fitting error and the SNR using the fitted sinc and the error column
    if 1 == 1:
        for mgo in range(len(newModels)):
            m = newModels[mgo]
            ff = newFeatures[mgo]
            #get the fitted parameters
            params = m.getFittedParameters()
            p0 = params[0]
            p1 = params[1]
            #get the error on the fitted parameters
            parametersStdDev = m.getFittedStdDev()
            measuredFreqErr = parametersStdDev[1]
            #create a sinc using the fitted parameters
            sinc = SincModel()
            sinc.setParameters(params)
            sinc = sinc(specIn.wave)
            #find the snr
            snrFull = sinc/specIn['error'].data
            pos = specIn.wave.where(ABS(sinc) == MAX(ABS(sinc)))
            snr = snrFull[pos]
            #get the continuum value at the same position
            contVal = continuum[pos]
            ###Remove any feature that isn't the right sign
            if sign == 'pos':
                #raise ValueError('sign:%s, snr:%s'%(sign,snr)) 
                if snr[0] < 0.:
                    fitter.removeModel(m)
                    continue
            if sign == 'neg':
                #raise ValueError('sign:%s, snr:%s'%(sign,snr))
                if snr[0] > 0.:
                    fitter.removeModel(m)
                    continue   
            ###Remove any feature that has snr < 500 CHANGE THIS TO 300
            if sign == 'neg':
                if snr[0] < -300:
                    fitter.removeModel(m)
            ###Remove any feature that has wandered well away
            if ABS(p1-ff) > checkValue:
                    fitter.removeModel(m)
            ###Remove any feature model with SNR < snrThreshold
#            if extraConstraints:
#                if ABS(snr) >= snrThreshold:
#                    ###Only keep features inside of the avoid regions
#                    if (p1 > (MIN(specIn.wave)+avoid)) & (p1 < (MAX(specIn.wave)-avoid)):
#                        ###Only keep features not wandered too far
#                        if ABS(p1-ff) < checkValue:
#                            newFeaturesResult.append(p1)
#                            newErrResult.append(measuredFreqErr)
#                            newSnrResult.append(snr)
#                            newPeakResult.append(p0)
#                            newModelResult.append(m)
#                            newContinuumResult.append(contVal)
#                        else:
#                            fitter.removeModel(m)
#                    else:
#                        fitter.removeModel(m)
#                else:
#                    fitter.removeModel(m)
#            else:
            if 1 == 1:
                #if ABS(snr) >= snrThreshold:
                ###Remove the snr check and leave that to the final SNR check at the end of the script
                if 1 == 1:
                    ###Only keep features inside of the avoid regions
                    ####Note that now all features will be inside of the avoid region
                    ####so really don't need this check
                    if (p1 > (MIN(specIn.wave)+avoid)) & (p1 < (MAX(specIn.wave)-avoid)):
                        ###Only keep features at a reasonable distance from the already accepted
                        ###This cuts down on doubles+ for partially resolved lines
                        if len(alreadyFoundFeatures) > 0:
                            min = MIN(ABS(p1-Double1d(alreadyFoundFeatures)))
                            pos = Double1d(alreadyFoundFeatures).where(ABS(p1-Double1d(alreadyFoundFeatures)) == MIN(ABS(p1-Double1d(alreadyFoundFeatures))))
                            #pos = Double1d(alreadyFoundFeatures).where(MIN(ABS(p1-Double1d(alreadyFoundFeatures))) <= 1.2)
                            #if pos.toInt1d().size > 0:
                            if ABS(p1-Double1d(alreadyFoundFeatures)[pos][0]) <=1.2:
                                #find the feature
                                moveSinc = sincModels[pos.toInt1d()[0]]
                                #set the position back to the one in alreadyFoundFeatures
                                moveParams = moveSinc.getFittedParameters()
                                moveParams[1] = alreadyFoundFeatures[pos.toInt1d()[0]]
                                moveSinc.setParameters(list(moveParams))
                                #And now try fixing the position
                                moveSinc.setFixed([1])
                                continue
                        newFeaturesResult.append(p1)
                        newErrResult.append(measuredFreqErr)
                        newSnrResult.append(snr)
                        newPeakResult.append(p0)
                        newModelResult.append(m)
                        newContinuumResult.append(contVal)
                        #now the model is accepted, limit its position
                        m.setLimits(1,p1-limitValue,p1+limitValue)
                    else:
                        fitter.removeModel(m)
                #else:
                #    fitter.removeModel(m)
    #Now repeat the fit with the vetted features
    fitter.doGlobalFit()
    fitter.residual()
    #Find the total model and residual
    totalModel = fitter.getTotalModel()
    residual = fitter.getResidual()
    #Update the fitted frequencies of already found features
    alreadyFoundFeatures = updateFeatureFrequencies(sincModels)
    #now update the flagging
    flags = flagSpec1d(specIn, alreadyFoundFeatures, flags.copy(), flagWidth, recreate = True)
    #Remove limits for kept features?
    #Outputs
    specOut = specIn.copy()
    specOut.flux = residual
    near = residualFlag2.wave.where((residualFlag2.wave < (MIN(residualFlag2.wave)+avoid)) | \
        (residualFlag2.wave > (MAX(residualFlag2.wave)-avoid)))
    residualFlag2.flag[near] = 1
    residualFlag3.flag[near] = 1
    #use flagWidth to flag the region round the features found
    if flagValue != None:
        for fq in newFeaturesResult:
            residualFlag2.flag[residualFlag2.wave.where(ABS(residualFlag2.wave - fq) <\
                flagWidth)] = flagValue
        #flags[specIn.meta['channelName'].value]['p1'].append(Double1d(newFeaturesResult))
        #flags[specIn.meta['channelName'].value]['fL'].append(Double1d(newFeaturesResult)-flagWidth)
        #flags[specIn.meta['channelName'].value]['fR'].append(Double1d(newFeaturesResult)+flagWidth)
        flags = flagSpec1d(specIn,newFeaturesResult,flags,flagWidth)
    return specOut, fitter, newFeaturesResult, newSnrResult, newPeakResult, newModelResult, \
           newErrResult, newContinuumResult, additMods,\
           flags, alreadyFoundFeatures, residualFlag2, residualFlag3
#           ,\
#           removeLimits

#main fuction to run the feature finding show
def catalogueFeatures(specIn, snrThresholds, signs, jumpThreshold, array, \
    flags, \
    mergeRadius=20.0,\
    resampleResolution=5.0, fwhm=None, flagWidths=None, polyOrder=3,\
    apod=False, snrCut=5.0, snrRange=[4.0,25.0], snrCap=500, subtractBaseline=True, \
    baselineOrder=3, checkCutResults=False, avoid=None, testGofFcut=False, \
    testCheckFit=False, testNeutralC=False, estimateRedshift=False,\
    testGofFcutLines=False, limitValue=2.0, maxIter=None, checkValue=2.0,\
    noExtraBox=False):
    """
    Iteratively finds features with SNRs above the respective snrThreshold.
    The final SNR are taken from the final fitted peaks and local RMS of the full 
    residual.
    Inputs:
        - specIn: spectrum for finding features
        - snrThresholds: successive SNR limits for finding features
        - signs: list specifying feature signs desired per SNR threshold
        - jumpThreshold: minimum size of jump to indicate findJumps should weight and flag a region [STDDEV]
        - mergeRadius: width within which features are merged by listFeatures [GHz]
        - resampleResolution:  resolution used by findJumps for resampling [GHz]
        - fwhm: FWHM of sinc profiles [GHz]. If None then actualResolution is used
        - flagWidths: list of half widths used by fitFeatures to flag features found [GHz]
        - polyOrder: degree of polynomial fitted to the continuum
        - apod: set to True if specIn is an apodized spectrum
        - snrCut: final SNR cut for features found  
        - snrRange: range either side of feature for final noise estimate [GHz]
        - snrCap: maximum SNR allowed when calculating the initial peak guess
        - subtractBaseline: set for an extra baseline subtraction when taking for the final SNR
        - baselineOrder: order of poly fitted for subtractBaseline
        - checkCutResults: set to a snrThreshold to stop fitting at that point
        - avoid: width of the band ends to exclude from the finding process [GHz]
        - testGofFcut: if set to True, GoF is run per snrThreshold
        - testCheckFit: set to True to check and update the fitted peak
        - testNeutralC: set to True to look for the CI feature
        - estimateReshift: set to True to estimate the redshift when looking for the CI feature
        - testGofFcutLines: testing GofF to remove features per SNR threshold
        - limitValue: limits the feature position during the fit
        - checkValue: if a feature has moved more than this, it is removed after the test fit
    Returns:
        - featuresFound: fitted frequencies of features found [GHz]
        - snrResult: final SNR of features (fitted peak/local noise)
        - snrIterResult: initial SNR of features, determined when found (peak/errorColumn)
        - cutFeaturesFound: fitted frequencies of features found > snrCut [GHz]
        - cutSnrResult: final SNR of features > snrCut (fitted peak/local noise)
        - cutSnrIterResult: initial SNR of features > snrCut, determined when found (peak/errorColumn)
        - snrThresholdsResult: snrThresholds at which features were found
        - cutSnrThresholdsResult: snrThresholds at which features > snrCut were found
        - finalTotalModel: final total fitted model (array)
        - finalResidual: final residual (array)
        - fitter: instance of SpectrumFitter used to fit all the models 
        - continuumModel: polynominal model fitted to the continuum
        - sincModels: final sinc models fitted to the features
        - cutSincModels: final sinc models fitted to the features > snrCut
        - featuresFoundErr: error on fitted feature frequencies
        - cutFeaturesFoundErr: error on fitted feature frequencies for features > snrCut
        - initialContinuumArray: initial polynominal fit to the continuum (array)
        - peaksAll: fitted peaks
        - cutPeaks: fitted peaks for features > snrCut
        - cutContinuumValues: corresponding continuum values to fitted peaks for features > snrCut
        - outputInitialFeatureFreq: initial fitted feature frequencies
        - cutOutputInitialFeatureFreq: initial fitted feature frequencies for features > snrCut
        - indexFeatures: feature number (used to create unique feature names per obsid)
        - cutIndexFeatures: feature number (used to create unique feature names per obsid, for features > snrCut)
        - threshSincs: sincs fitted at each snrThreshold cut (list of arrays)
        - threshResiduals: residual for each snrThreshold cut (list of arrays)
        - threshTotalModels: total model for each snrThreshold cut (list of arrays)    
        - allGofF: GoF for the features found per snrThreshold
        - allFeaturesForGofF: corresponding features for allGofF
        - allSnrForGofF: corresponding SNRs for allGofF
        - expectedFreq: expected frequency for a doppler shifted CO(7-6) or CO(10-9) feature
    """
    #If avoid is set, flag the ends of the input spectrum
#    if avoid:
#        near = specIn.wave.where((specIn.wave < (MIN(specIn.wave)+avoid)) | \
#        (specIn.wave > (MAX(specIn.wave)-avoid)))
#        specIn.flag[near] = 1
    #Preliminary detection of bright peaks, to flag before running the initial 
    #continuum fit *NOTE THAT NO ACTUAL FLAGGING TAKES PLACE*
    resampledSpec, weightedSpec = findJumps(specIn.copy(), jumpThreshold,\
                                  resampleResolution=resampleResolution, \
                                  weight=True, flag=False, equalWeight=apod)
    #Initial continuum fit
    #specOut is continuum subtracted and used in the first feature finder loop
    specOut, fitter, continuumModel, initialContinuumArray = removeContinuum(weightedSpec,\
                                   polyOrder)
    continuumArray = initialContinuumArray
    noiseFloor = MEAN(ABS(specOut.flux))*1.0
    #Try to estimate the redshift
    #if testNeutralC:
    #if estimateRedshift == 1:
    #    expectedFreq, neighbour, zEstimate = redshift1(specOut,array)
    #elif estimateRedshift == 2:
    #    expectedFreq, neighbour, zEstimate = redshift2(specOut,factor=20,mergeRadius=10)
    #else:
    #    expectedFreq, neighbour, zEstimate = -99, False, -99
    if checkCutResults:
        cutResiduals = Double1d()
        cutTotalModels = Double1d()
    #Number of snrThresholds to iterate over
    N = len(snrThresholds)
    if flagWidths == None:
        flagWidths = [0.0] * N
    #Set up the outputs
    threshSincs, threshResiduals ,threshTotalModels = [],[],[]
    featuresFound, outputInitialFeatureFreq, featuresFoundErr, snrResult = [],[],[],[]
    continuumPeakResult, continuumResult = [],[]
    additMods = []
    #Keep track of the SNR calculated during the iterations
    snrIterResult = []
    #Keep track of what iteration the feature was found in
    snrThresholdsResult = []
    #Keep track of the flagwidth used per feature
    flagWidthsAll  = []
    peaksAll ,sincModels ,cutSincModels = [],[],[]
    cutFeaturesFound ,cutFeaturesFoundErr = Double1d(), Double1d()
    cutSnrResult, cutSnrIterResult = Double1d(), Double1d()
    cutSnrResult2, cutSnrResult3 = Double1d(), Double1d()
    cutContinuumPeakValues ,cutContinuumValues, cutPeaks = Double1d(), Double1d(), Double1d()
    cutOutputInitialFeatureFreq, cutSnrThresholdsResult = Double1d(), Double1d()
    #Keep a list of the gOfF per threshold, along with the features and initial SNRs
    allGofF, allFeaturesForGofF, allSnrForGofF = [],[],[]
    ###Remove the fitting weights used for the initial continuum fit
    fitter.setWeight(None)
    #
    residualSpec2 = specIn.copy()
    residualSpec3 = specIn.copy()
    residualSpec2Out = specIn.copy()
    residualSpec3Out = specIn.copy()
    ###Search for features per snrThreshold
    for i in xrange(N):
        ###find the SNR using the error column
        snr = specOut.copy()
        snr.setFlux(snr.flux/snr['error'].data)
        snr.meta['obsid']=LongParameter(obsid)
        snr.meta['array']=StringParameter(array)
        ###find new features above the snrThreshold
        newFeatures, newSnr = listFeatures(snr, snrThresholds[i], flags, mergeRadius=mergeRadius, sign=signs[i])
        ###Fit the features *newFeatures (list of new features) is added to the existing features*
        #fitFeatures would need amending to fit only the new features
        #to the residual from the previous iteration.
        #try:
        specOut, fitter, newFeatures, newSnr, newPeaks, newModels, newFeaturesErr, \
            continuumValues, additMods,\
                         flags, featuresFound,\
                         residualSpec2out, residualSpec3out,\
                                    = fitFeatures(specIn,residualSpec2, residualSpec3, newFeatures, newSnr, continuumArray, \
                                      flags, sincModels, featuresFound, snrThresholdsResult, snrIterResult, \
                                      fwhm=fwhm, apod=apod, flagWidth=flagWidths[i], \
                                      snrThreshold=snrThresholds[i],fitter=fitter, snrCap=snrCap,\
                                      additMods=additMods,\
                                      noiseFloor=noiseFloor,avoid=avoid,limitValue=limitValue,\
                                      maxIter=maxIter,sign=signs[i])
        #except:
        #    continue
        #Update the continuum
        paramC = continuumModel.getFittedParameters()
        paramErrorC = continuumModel.getFittedStdDev()
        polyC = PolynomialModel(polyOrder)
        polyC.setParameters(paramC)
        continuumArray = polyC(specIn.wave)
        #Keep track of at which threshold the feature was found in
        thresholds = [snrThresholds[i]] * len(newFeatures)
        #Keep track of the flagWidths used for each feature
        widths = [flagWidths[i]] * len(newFeatures)
        #update flagging
        specIn.setFlag(specOut.flag)
        residualSpec2.setFlag(residualSpec2Out.flag)
        residualSpec3.setFlag(residualSpec3Out.flag)
        #get the total model and residual as arrays
        tempTotalMod = fitter.getTotalModel()
        tempResidual = fitter.getResidual()
        if checkCutResults:
            cutTotalModels.append(tempTotalMod)
            cutResiduals.append(tempResidual)
        #Add newFeatures and newPeaks to the previous lists
        featuresFound += newFeatures
        outputInitialFeatureFreq += newFeatures
        featuresFoundErr += newFeaturesErr        
        continuumPeakResult += (continuumValues + newPeaks)
        continuumResult +=continuumValues
        peaksAll += newPeaks
        snrIterResult += newSnr
        snrThresholdsResult += thresholds
        flagWidthsAll += widths
        sincModels += newModels
        threshResiduals.append(tempResidual)
        threshTotalModels.append(tempTotalMod)
        ###go through all the models and get the sincs
        sincResult = []
        for mm in sincModels:
           params = mm.getFittedParameters()
           sinc = SincModel()
           sinc.setParameters(params)
           sinc = sinc(specIn.wave)
           sincResult+=sinc
        threshSincs.append(sincResult)        
        #Get gOfF for each snrThreshold
        if testGofFcut:
            gOfFmidCut = gOfF(featuresFound, specIn, sincModels, rng=gOfFRng)
            #what about removing bad looking features for each cut?
            if testGofFcutLines:
                gOfFmidCut,featuresFound,featuresFoundErr,continuumPeakResult,\
                     continuumResult,peaksAll,snrIterResult,snrThresholdsResult,\
                     flagWidthsAll,sincModels
            allGofF.append(gOfFmidCut)
            allFeaturesForGofF.append(featuresFound)
            allSnrForGofF.append(snrIterResult)
        if checkCutResults == snrThresholds[i]:
            return featuresFound, snrResult, snrIterResult, cutTotalModels, \
                   cutResiduals, Double1d([0]), snrThresholdsResult, Double1d([0]), \
                   tempTotalMod, tempResidual, fitter, continuumModel, sincModels, \
                   Double1d([0]), featuresFoundErr, Double1d([0]), \
                   Double1d([0]), Double1d([0]), Double1d([0]), \
                   Double1d([0]), Double1d([0]), Double1d([0]), \
                   Double1d([0]), Double1d([0]),\
                   threshSincs,threshResiduals,threshTotalModels,\
                   allGofF, allFeaturesForGofF, allSnrForGofF,\
                   expectedFreq, zEstimate
        ###UPDATE FLAGWIDTHS
        #1. set flag to None
        #2. reflag avoid
        #3. use models to adjust the flagging per feature
#        checkFlag = specIn.flag.copy()
#        #specIn.setFlag(None)
#        specIn.setFlag(Int1d(specIn.flag.size))
#        if avoid:
#            near = specIn.wave.where((specIn.wave < (MIN(specIn.wave)+avoid)) | \
#            (specIn.wave > (MAX(specIn.wave)-avoid)))
#            specIn.flag[near] = 1
#        for gogo in range(len(featuresFound)):
#            modPara = sincModels[gogo].getFittedParameters()
#            feature = modPara[1]
#            specIn.flag[specIn.wave.where(ABS(specIn.wave - modPara[1]) < flagWidthsAll[gogo])] = 1
    #Get the total fitted model and the associated residual
    finalTotalModel = fitter.getTotalModel()
    finalResidual = fitter.getResidual()
    #Put the residual into a spectral container
    residualSpec = specIn.copy()
    residualSpec.flux = finalResidual
    residualSpec2.flux = finalResidual
    residualSpec3.flux = finalResidual
    #residualSpec2 = residualSpec.copy()
    ###Use the final fitted feature positions to update featuresFound
    #(outputInitialFeatureFreq gives the initial fitted feature position)
    #updateFeaturesFound = []
    #for gogo in range(len(featuresFound)):
    #    modPara = sincModels[gogo].getFittedParameters()
    #    updateFeaturesFound.append(modPara[1])
    #Now set an alternative .flag using the final features, for the final SNR estimate
    for gogo in range(len(featuresFound)):
        modPara = sincModels[gogo].getFittedParameters()
        feature = modPara[1]
        residualSpec3.flag[residualSpec3.wave.where(ABS(residualSpec3.wave - feature) < snrRange[0])] = 1
    #test = PlotXY()
    #test.addLayer(LayerXY(residualSpec.wave,residualSpec.flag))
    #stop
    #outputInitialFeatureFreq = featuresFound
    #featuresFound = updateFeaturesFound
    ###Index the features
    indexFeatures = String1d(range(len(featuresFound)))
    cutIndexFeatures = String1d()
    ###For each feature find the SNR using the full residual and fitted peaks
    for gogo in range(len(featuresFound)):
        feature = featuresFound[gogo]
        existingSnr = snrIterResult[gogo]
        if (feature > (MIN(specIn.wave)+avoid)) & (feature < (MAX(specIn.wave)-avoid)):
            featureInd = indexFeatures[gogo]
            err = featuresFoundErr[gogo]
            peak = peaksAll[gogo]
            snrIter = snrIterResult[gogo]
            threshold = snrThresholdsResult[gogo]
            continuumPeak = continuumPeakResult[gogo]
            continuumValue = continuumResult[gogo]
            singleMod = sincModels[gogo]
            initialFreq = outputInitialFeatureFreq[gogo]
            #find the final SNR using the fitted peak and local RMS of the full residal
            #no .flag
            snr = findResidualSnr(residualSpec, feature, peak, existingSnr, snrRange=snrRange, \
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Get the SNR using the old .flag
            snr2 = findResidualSnr(residualSpec2, feature, peak, existingSnr, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Get the SNR using the snrRange .flag
            snr3 = findResidualSnr(residualSpec3, feature, peak, existingSnr, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Keep the feature if it is above the SNR limit
            if ABS(snr3) >= snrCut:
                #for features with SNR > 10 remove features found that are 
                #within 1.5-2.1 GHz to the right (asymmetry wing fit removal)
                #The lesser feature should have << SNR compared to the main feature
                #print 'positive side'
                #****Need to update the model???????? and flagging??????
                if array[:3] == 'SLW':
                    pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) >= 1.5).\
                    and((feature-Double1d(featuresFound)) <= 2.1))
                elif array[:3] == 'SSW':
                    pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) > 0.0).\
                    and((feature-Double1d(featuresFound)) <= 2.1))
                if pos.length():
                    if len(pos.toInt1d()) > 1:
                        snrCheck = Double1d()
                        for fe in range(pos.toInt1d().size):
                            fefe = featuresFound[pos.toInt1d()[fe]]
                            pepe = peaksAll[pos.toInt1d()[fe]]
                            snr3temp = findResidualSnr(residualSpec3, fefe, pepe, 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                            snrCheck.append(snr3temp)
                        #find the strongest peak
                        pos = snrCheck.where(snrCheck == MAX(snrCheck))
                        maxSnr = snrCheck[pos]
                    else:
                        maxSnr = findResidualSnr(residualSpec3, Double1d(featuresFound)[pos][0], Double1d(peaksAll)[pos][0], 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)                            
                    if maxSnr > snr3:
                        if maxSnr > 10:
                            if (maxSnr/snr3 > 4):
                            #& (maxSnr/snr3 < 10):
                                continue
                ###Now similar for -ve features on the left of a strong line
                #for features with SNR > 10 remove features found that are 
                #within 1.5-2.1 GHz to the right (asymmetry wing fit removal)
                #The lesser feature should have << SNR compared to the main feature
                #print 'negative side'
                pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) < 0.0).\
                and((feature-Double1d(featuresFound)) >= -2.5))
                if pos.length():
                    if len(pos.toInt1d()) > 1:
                        snrCheck = Double1d()
                        for fe in range(pos.toInt1d().size):
                            fefe = featuresFound[pos.toInt1d()[fe]]
                            pepe = peaksAll[pos.toInt1d()[fe]]
                            snr3temp = findResidualSnr(residualSpec3, fefe, pepe, 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                            snrCheck.append(snr3temp)
                        #find the strongest peak
                        pos = snrCheck.where(ABS(snrCheck) == MAX(ABS(snrCheck)))
                        maxSnr = snrCheck[pos]
                    else:
                        maxSnr = findResidualSnr(residualSpec3, Double1d(featuresFound)[pos][0], Double1d(peaksAll)[pos][0], 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                    if maxSnr > ABS(snr3):
                        if maxSnr > 10:
                            if (ABS(maxSnr/snr3) > 4):
                                #& (maxSnr/snr3 < 10):
                                continue  
                cutFeaturesFound.append(feature)
                cutIndexFeatures.append(featureInd)
                cutFeaturesFoundErr.append(err)
                cutSnrResult.append(snr)
                cutSnrResult2.append(snr2)
                cutSnrResult3.append(snr3)
                cutSnrIterResult.append(snrIter)
                cutSnrThresholdsResult.append(threshold)
                cutPeaks.append(peak)
                cutContinuumPeakValues.append(continuumPeak)
                cutContinuumValues.append(continuumValue)
                cutSincModels.append(singleMod)
                cutOutputInitialFeatureFreq.append(initialFreq)
            snrResult += [snr]
    return featuresFound, snrResult, snrIterResult, cutFeaturesFound, cutSnrResult, \
             cutSnrIterResult, snrThresholdsResult, cutSnrThresholdsResult, \
             finalTotalModel, finalResidual, fitter, continuumModel, sincModels, \
             cutSincModels, featuresFoundErr, cutFeaturesFoundErr,\
             initialContinuumArray, peaksAll, cutPeaks, \
             cutContinuumValues, outputInitialFeatureFreq, cutOutputInitialFeatureFreq, \
             indexFeatures, cutIndexFeatures,\
             threshSincs,threshResiduals,threshTotalModels,\
             allGofF, allFeaturesForGofF, allSnrForGofF,\
             flags, cutSnrResult2, cutSnrResult3,\
             residualSpec,residualSpec2,residualSpec3

def catalogueFeaturesWithPp(specIn, snrThresholds, signs, jumpThreshold, array, \
    flags, \
    mergeRadius=20.0,\
    resampleResolution=5.0, fwhm=None, flagWidths=None, polyOrder=3,\
    apod=False, snrCut=5.0, snrRange=[4.0,25.0], snrCap=500, subtractBaseline=True, \
    baselineOrder=3, checkCutResults=False, avoid=None, testGofFcut=False, \
    testCheckFit=False, testNeutralC=False, estimateRedshift=False,\
    testGofFcutLines=False, limitValue=2.0, maxIter=None, checkValue=2.0,\
    noExtraBox=False, ppTable=False):
    """
    Iteratively finds features with SNRs above the respective snrThreshold.
    The final SNR are taken from the final fitted peaks and local RMS of the full 
    residual.
    Inputs:
        - specIn: spectrum for finding features
        - snrThresholds: successive SNR limits for finding features
        - signs: list specifying feature signs desired per SNR threshold
        - jumpThreshold: minimum size of jump to indicate findJumps should weight and flag a region [STDDEV]
        - mergeRadius: width within which features are merged by listFeatures [GHz]
        - resampleResolution:  resolution used by findJumps for resampling [GHz]
        - fwhm: FWHM of sinc profiles [GHz]. If None then actualResolution is used
        - flagWidths: list of half widths used by fitFeatures to flag features found [GHz]
        - polyOrder: degree of polynomial fitted to the continuum
        - apod: set to True if specIn is an apodized spectrum
        - snrCut: final SNR cut for features found  
        - snrRange: range either side of feature for final noise estimate [GHz]
        - snrCap: maximum SNR allowed when calculating the initial peak guess
        - subtractBaseline: set for an extra baseline subtraction when taking for the final SNR
        - baselineOrder: order of poly fitted for subtractBaseline
        - checkCutResults: set to a snrThreshold to stop fitting at that point
        - avoid: width of the band ends to exclude from the finding process [GHz]
        - testGofFcut: if set to True, GoF is run per snrThreshold
        - testCheckFit: set to True to check and update the fitted peak
        - testNeutralC: set to True to look for the CI feature
        - estimateReshift: set to True to estimate the redshift when looking for the CI feature
        - testGofFcutLines: testing GofF to remove features per SNR threshold
        - limitValue: limits the feature position during the fit
        - checkValue: if a feature has moved more than this, it is removed after the test fit
    Returns:
        - featuresFound: fitted frequencies of features found [GHz]
        - snrResult: final SNR of features (fitted peak/local noise)
        - snrIterResult: initial SNR of features, determined when found (peak/errorColumn)
        - cutFeaturesFound: fitted frequencies of features found > snrCut [GHz]
        - cutSnrResult: final SNR of features > snrCut (fitted peak/local noise)
        - cutSnrIterResult: initial SNR of features > snrCut, determined when found (peak/errorColumn)
        - snrThresholdsResult: snrThresholds at which features were found
        - cutSnrThresholdsResult: snrThresholds at which features > snrCut were found
        - finalTotalModel: final total fitted model (array)
        - finalResidual: final residual (array)
        - fitter: instance of SpectrumFitter used to fit all the models 
        - continuumModel: polynominal model fitted to the continuum
        - sincModels: final sinc models fitted to the features
        - cutSincModels: final sinc models fitted to the features > snrCut
        - featuresFoundErr: error on fitted feature frequencies
        - cutFeaturesFoundErr: error on fitted feature frequencies for features > snrCut
        - initialContinuumArray: initial polynominal fit to the continuum (array)
        - peaksAll: fitted peaks
        - cutPeaks: fitted peaks for features > snrCut
        - cutContinuumValues: corresponding continuum values to fitted peaks for features > snrCut
        - outputInitialFeatureFreq: initial fitted feature frequencies
        - cutOutputInitialFeatureFreq: initial fitted feature frequencies for features > snrCut
        - indexFeatures: feature number (used to create unique feature names per obsid)
        - cutIndexFeatures: feature number (used to create unique feature names per obsid, for features > snrCut)
        - threshSincs: sincs fitted at each snrThreshold cut (list of arrays)
        - threshResiduals: residual for each snrThreshold cut (list of arrays)
        - threshTotalModels: total model for each snrThreshold cut (list of arrays)    
        - allGofF: GoF for the features found per snrThreshold
        - allFeaturesForGofF: corresponding features for allGofF
        - allSnrForGofF: corresponding SNRs for allGofF
        - expectedFreq: expected frequency for a doppler shifted CO(7-6) or CO(10-9) feature
    """
    #Preliminary detection of bright peaks, to flag before running the initial 
    #continuum fit *NOTE THAT NO ACTUAL FLAGGING TAKES PLACE*
    resampledSpec, weightedSpec = findJumps(specIn.copy(), jumpThreshold,\
                                  resampleResolution=resampleResolution, \
                                  weight=True, flag=False, equalWeight=apod)
    #Initial continuum fit
    #specOut is continuum subtracted and used in the first feature finder loop
    specOut, fitter, continuumModel, initialContinuumArray = removeContinuum(weightedSpec,\
                                   polyOrder)
    continuumArray = initialContinuumArray
    noiseFloor = MEAN(ABS(specOut.flux))*1.0
    if checkCutResults:
        cutResiduals = Double1d()
        cutTotalModels = Double1d()
    #Number of snrThresholds to iterate over
    N = len(snrThresholds)
    if flagWidths == None:
        flagWidths = [0.0] * N
    #Set up the outputs
    threshSincs, threshResiduals ,threshTotalModels = [],[],[]
    featuresFound, outputInitialFeatureFreq, featuresFoundErr, snrResult = [],[],[],[]
    continuumPeakResult, continuumResult = [],[]
    additMods = []
    #Keep track of the SNR calculated during the iterations
    snrIterResult = []
    #Keep track of what iteration the feature was found in
    snrThresholdsResult = []
    #Keep track of the flagwidth used per feature
    flagWidthsAll  = []
    peaksAll ,sincModels ,cutSincModels = [],[],[]
    cutFeaturesFound ,cutFeaturesFoundErr = Double1d(), Double1d()
    cutSnrResult, cutSnrIterResult = Double1d(), Double1d()
    cutSnrResult2, cutSnrResult3 = Double1d(), Double1d()
    cutContinuumPeakValues ,cutContinuumValues, cutPeaks = Double1d(), Double1d(), Double1d()
    cutOutputInitialFeatureFreq, cutSnrThresholdsResult = Double1d(), Double1d()
    #Keep a list of the gOfF per threshold, along with the features and initial SNRs
    allGofF, allFeaturesForGofF, allSnrForGofF = [],[],[]
    ###Remove the fitting weights used for the initial continuum fit
    fitter.setWeight(None)
    #
    residualSpec2 = specIn.copy()
    residualSpec3 = specIn.copy()
    residualSpec2Out = specIn.copy()
    residualSpec3Out = specIn.copy()
    ###Search for features per snrThreshold
    if ppTable:
#        ppFreq = ppTable['data']['frequency'].data.copy()
#        ppSnr = ppTable['data']['SNR'].data.copy()
#        ppDet = ppTable['data']['detector'].data.copy()
#        ppPeak = ppTable['data']['fluxSub'].data.copy()
        ppFreq = ppTable['frequency'].data.copy()
        ppSnr = ppTable['SNR'].data.copy()
        ppDet = ppTable['detector'].data.copy()
        ppPeak = ppTable['fluxSub'].data.copy()
        posPpDet = ppDet.where(ppDet == array[:3])
        #print posPpDet, array
        if posPpDet.toInt1d().size > 0:
            #print 'Using PP results'
            ppFreq = ppFreq[posPpDet]
            ppSnr = ppSnr[posPpDet]*(5/3.)
            ppPeak = ppPeak[posPpDet]
    #
    for i in xrange(N):
        ###find the SNR using the error column
        snr = specOut.copy()
        snr.setFlux(snr.flux/snr['error'].data)
        snr.meta['obsid']=LongParameter(obsid)
        snr.meta['array']=StringParameter(array)
        ###find new features above the snrThreshold
        if not ppTable:
            newFeatures, newSnr = listFeatures(snr, snrThresholds[i], flags, mergeRadius=mergeRadius, sign=signs[i])
        else:
            posPp = ppSnr.where(ppSnr > snrThresholds[i])
            if posPp.toInt1d().size > 0:
                newFeatures = list(ppFreq[posPp].copy())
                newSnr = list(ppPeak[posPp].copy())
                ppSnr[posPp] = 0
            else:
                newFeatures = []
                newSnr = []
                #print newFeatures,newSnr
        ###Fit the features *newFeatures (list of new features) is added to the existing features*
        #fitFeatures would need amending to fit only the new features
        #to the residual from the previous iteration.
        #try:
        specOut, fitter, newFeatures, newSnr, newPeaks, newModels, newFeaturesErr, \
            continuumValues, additMods,\
                         flags, featuresFound,\
                         residualSpec2out, residualSpec3out,\
                                    = fitFeatures(specIn,residualSpec2, residualSpec3, newFeatures, newSnr, continuumArray, \
                                      flags, sincModels, featuresFound, snrThresholdsResult, snrIterResult, \
                                      fwhm=fwhm, apod=apod, flagWidth=flagWidths[i], \
                                      snrThreshold=snrThresholds[i],fitter=fitter, snrCap=snrCap,\
                                      additMods=additMods,\
                                      noiseFloor=noiseFloor,avoid=avoid,limitValue=limitValue,\
                                      maxIter=maxIter,sign=signs[i],ppTable=ppTable)
        #except:
        #    continue
        #Update the continuum
        paramC = continuumModel.getFittedParameters()
        paramErrorC = continuumModel.getFittedStdDev()
        polyC = PolynomialModel(polyOrder)
        polyC.setParameters(paramC)
        continuumArray = polyC(specIn.wave)
        #Keep track of at which threshold the feature was found in
        thresholds = [snrThresholds[i]] * len(newFeatures)
        #Keep track of the flagWidths used for each feature
        widths = [flagWidths[i]] * len(newFeatures)
        #update flagging
        specIn.setFlag(specOut.flag)
        residualSpec2.setFlag(residualSpec2Out.flag)
        residualSpec3.setFlag(residualSpec3Out.flag)
        #get the total model and residual as arrays
        tempTotalMod = fitter.getTotalModel()
        tempResidual = fitter.getResidual()
        if checkCutResults:
            cutTotalModels.append(tempTotalMod)
            cutResiduals.append(tempResidual)
        #Add newFeatures and newPeaks to the previous lists
        featuresFound += newFeatures
        outputInitialFeatureFreq += newFeatures
        featuresFoundErr += newFeaturesErr        
        continuumPeakResult += (continuumValues + newPeaks)
        continuumResult +=continuumValues
        peaksAll += newPeaks
        snrIterResult += newSnr
        snrThresholdsResult += thresholds
        flagWidthsAll += widths
        sincModels += newModels
        threshResiduals.append(tempResidual)
        threshTotalModels.append(tempTotalMod)
        ###go through all the models and get the sincs
        sincResult = []
        for mm in sincModels:
           params = mm.getFittedParameters()
           sinc = SincModel()
           sinc.setParameters(params)
           sinc = sinc(specIn.wave)
           sincResult+=sinc
        threshSincs.append(sincResult)        
        #Get gOfF for each snrThreshold
        if testGofFcut:
            gOfFmidCut = gOfF(featuresFound, specIn, sincModels, rng=gOfFRng)
            #what about removing bad looking features for each cut?
            if testGofFcutLines:
                gOfFmidCut,featuresFound,featuresFoundErr,continuumPeakResult,\
                     continuumResult,peaksAll,snrIterResult,snrThresholdsResult,\
                     flagWidthsAll,sincModels
            allGofF.append(gOfFmidCut)
            allFeaturesForGofF.append(featuresFound)
            allSnrForGofF.append(snrIterResult)
        if checkCutResults == snrThresholds[i]:
            return featuresFound, snrResult, snrIterResult, cutTotalModels, \
                   cutResiduals, Double1d([0]), snrThresholdsResult, Double1d([0]), \
                   tempTotalMod, tempResidual, fitter, continuumModel, sincModels, \
                   Double1d([0]), featuresFoundErr, Double1d([0]), \
                   Double1d([0]), Double1d([0]), Double1d([0]), \
                   Double1d([0]), Double1d([0]), Double1d([0]), \
                   Double1d([0]), Double1d([0]),\
                   threshSincs,threshResiduals,threshTotalModels,\
                   allGofF, allFeaturesForGofF, allSnrForGofF,\
                   expectedFreq, zEstimate
        ###UPDATE FLAGWIDTHS
        #1. set flag to None
        #2. reflag avoid
        #3. use models to adjust the flagging per feature
#        checkFlag = specIn.flag.copy()
#        #specIn.setFlag(None)
#        specIn.setFlag(Int1d(specIn.flag.size))
#        if avoid:
#            near = specIn.wave.where((specIn.wave < (MIN(specIn.wave)+avoid)) | \
#            (specIn.wave > (MAX(specIn.wave)-avoid)))
#            specIn.flag[near] = 1
#        for gogo in range(len(featuresFound)):
#            modPara = sincModels[gogo].getFittedParameters()
#            feature = modPara[1]
#            specIn.flag[specIn.wave.where(ABS(specIn.wave - modPara[1]) < flagWidthsAll[gogo])] = 1
    #Get the total fitted model and the associated residual
    finalTotalModel = fitter.getTotalModel()
    finalResidual = fitter.getResidual()
    #Put the residual into a spectral container
    residualSpec = specIn.copy()
    residualSpec.flux = finalResidual
    residualSpec2.flux = finalResidual
    residualSpec3.flux = finalResidual
    #residualSpec2 = residualSpec.copy()
    ###Use the final fitted feature positions to update featuresFound
    #(outputInitialFeatureFreq gives the initial fitted feature position)
    #updateFeaturesFound = []
    #for gogo in range(len(featuresFound)):
    #    modPara = sincModels[gogo].getFittedParameters()
    #    updateFeaturesFound.append(modPara[1])
    #Now set an alternative .flag using the final features, for the final SNR estimate
    for gogo in range(len(featuresFound)):
        modPara = sincModels[gogo].getFittedParameters()
        feature = modPara[1]
        residualSpec3.flag[residualSpec3.wave.where(ABS(residualSpec3.wave - feature) < snrRange[0])] = 1
    #test = PlotXY()
    #test.addLayer(LayerXY(residualSpec.wave,residualSpec.flag))
    #stop
    #outputInitialFeatureFreq = featuresFound
    #featuresFound = updateFeaturesFound
    ###Index the features
    indexFeatures = String1d(range(len(featuresFound)))
    cutIndexFeatures = String1d()
    ###For each feature find the SNR using the full residual and fitted peaks
    for gogo in range(len(featuresFound)):
        feature = featuresFound[gogo]
        existingSnr = snrIterResult[gogo]
        if (feature > (MIN(specIn.wave)+avoid)) & (feature < (MAX(specIn.wave)-avoid)):
            featureInd = indexFeatures[gogo]
            err = featuresFoundErr[gogo]
            peak = peaksAll[gogo]
            snrIter = snrIterResult[gogo]
            threshold = snrThresholdsResult[gogo]
            continuumPeak = continuumPeakResult[gogo]
            continuumValue = continuumResult[gogo]
            singleMod = sincModels[gogo]
            initialFreq = outputInitialFeatureFreq[gogo]
            #find the final SNR using the fitted peak and local RMS of the full residal
            #no .flag
            snr = findResidualSnr(residualSpec, feature, peak, existingSnr, snrRange=snrRange, \
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Get the SNR using the old .flag
            snr2 = findResidualSnr(residualSpec2, feature, peak, existingSnr, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Get the SNR using the snrRange .flag
            snr3 = findResidualSnr(residualSpec3, feature, peak, existingSnr, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
            #Keep the feature if it is above the SNR limit
            #print feature,existingSnr,snr3
            if ABS(snr3) >= snrCut:
                #for features with SNR > 10 remove features found that are 
                #within 1.5-2.1 GHz to the right (asymmetry wing fit removal)
                #The lesser feature should have << SNR compared to the main feature
                #print 'positive side'
                #****Need to update the model???????? and flagging??????
                if array[:3] == 'SLW':
                    pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) >= 1.5).\
                    and((feature-Double1d(featuresFound)) <= 2.1))
                elif array[:3] == 'SSW':
                    pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) > 0.0).\
                    and((feature-Double1d(featuresFound)) <= 2.1))
                if pos.length():
                    if len(pos.toInt1d()) > 1:
                        snrCheck = Double1d()
                        for fe in range(pos.toInt1d().size):
                            fefe = featuresFound[pos.toInt1d()[fe]]
                            pepe = peaksAll[pos.toInt1d()[fe]]
                            snr3temp = findResidualSnr(residualSpec3, fefe, pepe, 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                            snrCheck.append(snr3temp)
                        #find the strongest peak
                        pos = snrCheck.where(snrCheck == MAX(snrCheck))
                        maxSnr = snrCheck[pos]
                    else:
                        maxSnr = findResidualSnr(residualSpec3, Double1d(featuresFound)[pos][0], Double1d(peaksAll)[pos][0], 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)                            
                    if maxSnr > snr3:
                        if maxSnr > 10:
                            if (maxSnr/snr3 > 4):
                            #& (maxSnr/snr3 < 10):
                                continue
                ###Now similar for -ve features on the left of a strong line
                #for features with SNR > 10 remove features found that are 
                #within 1.5-2.1 GHz to the right (asymmetry wing fit removal)
                #The lesser feature should have << SNR compared to the main feature
                #print 'negative side'
                pos = Double1d(featuresFound).where(((feature-Double1d(featuresFound)) < 0.0).\
                and((feature-Double1d(featuresFound)) >= -2.5))
                if pos.length():
                    if len(pos.toInt1d()) > 1:
                        snrCheck = Double1d()
                        for fe in range(pos.toInt1d().size):
                            fefe = featuresFound[pos.toInt1d()[fe]]
                            pepe = peaksAll[pos.toInt1d()[fe]]
                            snr3temp = findResidualSnr(residualSpec3, fefe, pepe, 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                            snrCheck.append(snr3temp)
                        #find the strongest peak
                        pos = snrCheck.where(ABS(snrCheck) == MAX(ABS(snrCheck)))
                        maxSnr = snrCheck[pos]
                    else:
                        maxSnr = findResidualSnr(residualSpec3, Double1d(featuresFound)[pos][0], Double1d(peaksAll)[pos][0], 0.0, snrRange = snrRange,\
                                 ignoreFlag=1, subtractBaseline=True, \
                                 baselineOrder=baselineOrder, noExtraBox=noExtraBox)
                    if maxSnr > ABS(snr3):
                        if maxSnr > 10:
                            if (ABS(maxSnr/snr3) > 4):
                                #& (maxSnr/snr3 < 10):
                                continue  
                cutFeaturesFound.append(feature)
                cutIndexFeatures.append(featureInd)
                cutFeaturesFoundErr.append(err)
                cutSnrResult.append(snr)
                cutSnrResult2.append(snr2)
                cutSnrResult3.append(snr3)
                cutSnrIterResult.append(snrIter)
                cutSnrThresholdsResult.append(threshold)
                cutPeaks.append(peak)
                cutContinuumPeakValues.append(continuumPeak)
                cutContinuumValues.append(continuumValue)
                cutSincModels.append(singleMod)
                cutOutputInitialFeatureFreq.append(initialFreq)
            snrResult += [snr]
    return featuresFound, snrResult, snrIterResult, cutFeaturesFound, cutSnrResult, \
             cutSnrIterResult, snrThresholdsResult, cutSnrThresholdsResult, \
             finalTotalModel, finalResidual, fitter, continuumModel, sincModels, \
             cutSincModels, featuresFoundErr, cutFeaturesFoundErr,\
             initialContinuumArray, peaksAll, cutPeaks, \
             cutContinuumValues, outputInitialFeatureFreq, cutOutputInitialFeatureFreq, \
             indexFeatures, cutIndexFeatures,\
             threshSincs,threshResiduals,threshTotalModels,\
             allGofF, allFeaturesForGofF, allSnrForGofF,\
             flags, cutSnrResult2, cutSnrResult3,\
             residualSpec,residualSpec2,residualSpec3
