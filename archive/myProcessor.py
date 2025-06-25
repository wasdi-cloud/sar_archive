import wasdi
from datetime import datetime
from datetime import timedelta
import numpy

def deletePermanentWaterTiles(dLatN, dLonW, dLatS, dLonE, dLatStep, dLonStep, sMosaicFileNameBase):

    wasdi.wasdiLog("Deleting Permanent Water Tiles")

    iLatStepCount = 0
    iLonStepCount = 0

    for dLon in numpy.arange(dLonW, dLonE, dLonStep):
        for dLat in numpy.arange(dLatS, dLatN, dLatStep):
            sTileFileName = sMosaicFileNameBase + 'PERMWATER_PW__' + str(iLatStepCount) + '_' + str(iLonStepCount) + '.tif'
            wasdi.deleteProduct(sTileFileName)
            iLatStepCount = iLatStepCount + 1
        iLatStepCount = 0
        iLonStepCount = iLonStepCount + 1

def findInvolvedOrbits(sStartDate, sEndDate, oBBox, sPROVIDER, iMaxSearchDays, sPlatformNameFilter):
    # Start Date
    oStartDate = datetime.today()

    # Cast dates from params
    try:
        oStartDate = datetime.strptime(sStartDate, '%Y-%m-%d')
    except:
        wasdi.wasdiLog('Start Date not valid, assuming today')

    # End Date
    oEndDate = datetime.today()

    # Cast dates from params
    try:
        oEndDate = datetime.strptime(sEndDate, '%Y-%m-%d')
    except:
        wasdi.wasdiLog('End Date not valid, assuming today')

    oDelta = oEndDate - oStartDate

    oOriginalStartDate = oStartDate

    if oDelta.days>iMaxSearchDays:
        oTimeDeltaBack = timedelta(days=iMaxSearchDays)
        oStartDate = oEndDate-oTimeDeltaBack
        sStartDate = oStartDate.strftime("%Y-%m-%d")

    wasdi.wasdiLog("Search orbits Period = " + sStartDate + " - " + sEndDate)
    #aoReturnList = wasdi.searchEOImages("S1", sStartDate, sEndDate, sProductType="GRD", oBoundingBox=oBBox,
    #                                    sProvider=sPROVIDER)

    # if needed, filter based on platform name(s) for the search of all the images, for all the orbits
    if sPlatformNameFilter == "":
        wasdi.wasdiLog("No filter applied to the search of S1 images")
        aoReturnList = wasdi.searchEOImages("S1", sStartDate, sEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

    else:
        wasdi.wasdiLog("Only searching for S1 images of platforms: " + str(sPlatformNameFilter))
        #we need to remove whitespaces!!!!!!!!!!!!!
        asPlatformFilter = [x.strip() for x in sPlatformNameFilter.split(',')]
        aoReturnList = []
        aoReturnListAllPlatforms = wasdi.searchEOImages("S1", sStartDate, sEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

        for oReturnListAllPlatforms in aoReturnListAllPlatforms:
            if oReturnListAllPlatforms["fileName"].startswith(tuple(asPlatformFilter)):
                aoReturnList.append(oReturnListAllPlatforms)


    if oDelta.days>iMaxSearchDays:
        oTimeDeltaBack = timedelta(days=iMaxSearchDays)
        oNewEndDate = oOriginalStartDate+oTimeDeltaBack
        sStartDate = oOriginalStartDate.strftime("%Y-%m-%d")
        sNewEndDate = oNewEndDate.strftime("%Y-%m-%d")

        wasdi.wasdiLog("The period of the archive is very long, we search orbits also from start Period = " +sStartDate + " - " + sNewEndDate)

        #aoReturnList2 = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

        # if needed, filter based on platform name(s) for the search of all the images, for all the orbits
        if sPlatformNameFilter == "":
            wasdi.wasdiLog("No filter applied to the search of S1 images")
            aoReturnList2 = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

        else:
            wasdi.wasdiLog("Only searching for S1 images of platforms: " + str(sPlatformNameFilter))
            # we need to remove whitespaces!!!!!!!!!!!!!
            asPlatformFilter = [x.strip() for x in sPlatformNameFilter.split(',')]
            aoReturnList2 = []
            aoReturnList2AllPlatforms = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

            for oReturnList2AllPlatforms in aoReturnList2AllPlatforms:
                if oReturnList2AllPlatforms["fileName"].startswith(tuple(asPlatformFilter)):
                    aoReturnList2.append(oReturnList2AllPlatforms)

        if aoReturnList2 is not None:
            if len(aoReturnList2) > 0:
                aoReturnList.extend(aoReturnList2)

    if oEndDate.year>=2016:
        sStartDate = "2016-09-01"
        sNewEndDate = "2017-01-01"

        wasdi.wasdiLog("The period includes 2016 or later: we force a search also to verify if S1B is an alternative. Period = " +sStartDate + " - " + sNewEndDate)
        #aoReturnList3 = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

        # if needed, filter based on platform name(s) for the search of all the images, for all the orbits
        if sPlatformNameFilter == "":
            wasdi.wasdiLog("No filter applied to the search of S1 images")
            aoReturnList3 = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

        else:
            wasdi.wasdiLog("Only searching for S1 images of platforms: " + str(sPlatformNameFilter))
            # we need to remove whitespaces!!!!!!!!!!!!!
            asPlatformFilter = [x.strip() for x in sPlatformNameFilter.split(',')]
            aoReturnList3 = []
            aoReturnList3AllPlatforms = wasdi.searchEOImages("S1", sStartDate, sNewEndDate, sProductType="GRD", oBoundingBox=oBBox, sProvider=sPROVIDER)

            for oReturnList3AllPlatforms in aoReturnList3AllPlatforms:
                if oReturnList3AllPlatforms["fileName"].startswith(tuple(asPlatformFilter)):
                    aoReturnList3.append(oReturnList3AllPlatforms)

        if aoReturnList3 is not None:
            if len(aoReturnList3) > 0:
                aoReturnList.extend(aoReturnList3)

    # Find the list of unique orbits
    aiOrbits = []

    if len(aoReturnList) > 0:
        for oFoundImage in aoReturnList:
            iOrbit = oFoundImage["properties"]["relativeorbitnumber"]

            if iOrbit not in aiOrbits:
                aiOrbits.append(iOrbit)

    sOrbits = ""
    for iOrbit in aiOrbits:
        sOrbits += iOrbit + ","

    sOrbits = sOrbits[0:len(sOrbits) - 1]

    wasdi.wasdiLog("Automatic Orbits Found: [" + sOrbits + "]")

    return sOrbits

def run():
    wasdi.wasdiLog('HASARD Archive Generator v.3.4.3')

    sArchiveStartDate = wasdi.getParameter('ARCHIVE_START_DATE', "2019-05-01")
    sArchiveEndDate = wasdi.getParameter('ARCHIVE_END_DATE', "2019-05-03")

    # Read Bounding Box
    oBbox = wasdi.getParameter('BBOX')
    # Read Orbits
    sOrbits = wasdi.getParameter('ORBITS', "")
    # Read Search range in days
    iLastDays = wasdi.getParameter('LASTDAYS', 0)
    # Read Base Mosaic Name
    sMosaicFileNameBase = wasdi.getParameter('MOSAICBASENAME', 'mosaic')
    # Read suffix to apply to the daily flood maps
    #sSuffix = wasdi.getParameter('SUFFIX', '_flood.tif')
    sSuffix = wasdi.getParameter('SUFFIX', 'flood.tif')

    if sSuffix is None:
        sSuffix = ""

    if "_" in sSuffix:
        sSuffix = sSuffix.replace("_","")

    if sSuffix == "":
        wasdi.wasdiLog("Setting default suffix flood.tif")
        sSuffix = "flood.tif"
    elif not sSuffix.endswith(".tif"):
        sSuffix += ".tif"

    # Grid Lat, Lon Step
    sGridStep = wasdi.getParameter('GRIDSTEP', '1,1')
    # Delete flag
    bDelete = wasdi.getParameter('DELETE', False)
    # No Data Value
    iNODATAVALUE = wasdi.getParameter('NODATAVALUE', -9999)
    # Ignore Value
    iINPUTIGNOREVALUE = wasdi.getParameter('INPUTIGNOREVALUE', 0)
    # Mosaic No Data Value
    iMOSAICNODATAVALUE = wasdi.getParameter('MOSAICNODATAVALUE', 255)
    # Mosaic Input Ignore Value
    iMOSAICINPUTIGNOREVALUE = wasdi.getParameter('MOSAICINPUTIGNOREVALUE', 255)
    # Max Reference Tile Age in Days
    iMaxReferenceTileAge = wasdi.getParameter('MAXREFERENCETILEAGE', 90)
    # Data Provider
    sPROVIDER = wasdi.getParameter('PROVIDER', 'AUTO')
    # Flag to force re-run
    bForceReRun = wasdi.getParameter('FORCE_RE_RUN', False)
    # Read Simulate Flag
    bSimulate = wasdi.getParameter('SIMULATE', False)
    # Shape file mask
    sShapeMaskFile = wasdi.getParameter('SHAPEFILEMASK')
    # Read the workflow for preprocessing
    sWorkflow = wasdi.getParameter('PREPROCWORKFLOW', 'LISTSinglePreproc2')
    # Copernicus Dem Water Map File
    sCopDemMap = wasdi.getParameter("copdem_wm", None)
    # Apply the 3 state map conversion or not?
    bApplyMapConversion = wasdi.getParameter("APPLYMAPCONVERSION", False)
    # Minimum % intersection of a Sentinel image to be considered
    iMinimumPercentage = wasdi.getParameter("MINCOVERAGE", 0)
    # Flag to know if we need to move the result to sftp server or not
    bMoveToSftp = wasdi.getParameter("MOVETOSFTP", False)
    # Max number of days to use in the search to find automatic orbits
    iMaxSearchDays = wasdi.getParameter("MAX_SEARCH_DAYS_FOR_ORBITS", 120)
    # Flag to know if we need to apply a permanent water map
    bApplyWaterMap = wasdi.getParameter('APPLY_WATER_MAP', False)
    # Name of the Permanent Water Map Name. If none, it will take the Copernicus one
    sPermanentWaterMapName = wasdi.getParameter('PERMANENT_WATER_MAP_NAME', "")
    # Flood Detection
    sFloodDetectionProcessor = wasdi.getParameter("DETECTION_PROCESSOR", "HASARD")
    # LU LC Map to use in case we use AUTOWADE
    sLULCMap = wasdi.getParameter("LULC_IMAGE", None)
    # AUTOWADE Min Cluster
    iMinCluster = wasdi.getParameter("MIN_CLUST_N", 6)
    # AUTOWADE Filter Size
    iFilterSize = wasdi.getParameter("FILTER_SIZE", 3)
    # Filter on Platform names: S1A, S1B, S2C
    sPlatformFilter = wasdi.getParameter("PLATFORM_FILTER", "")

    if sFloodDetectionProcessor not in ["HASARD", "AUTOWADE"]:
        sFloodDetectionProcessor = "HASARD"

    bUseHasard = True

    if sFloodDetectionProcessor == "AUTOWADE":
        wasdi.wasdiLog("Using Autowade as Detection Processor")
        bUseHasard = False
    else:
        wasdi.wasdiLog("Using HASARD as Detection Processor")

    wasdi.wasdiLog('Start SAR Archive Generation from ' + sArchiveStartDate + ' to ' + sArchiveEndDate)

    if "_" in sMosaicFileNameBase:
        wasdi.wasdiLog('Remove _ from base name ')
        sMosaicFileNameBase = sMosaicFileNameBase.replace("_", "")

    if " " in sMosaicFileNameBase:
        wasdi.wasdiLog('Remove spaces from base name ')
        sMosaicFileNameBase = sMosaicFileNameBase.replace(" ", "")

    oStartDay = datetime.today()
    oEndDay = datetime.today()

    try:
        oStartDay = datetime.strptime(sArchiveStartDate, '%Y-%m-%d')
    except:
        wasdi.wasdiLog('Start Date not valid, assuming today')

    try:
        oEndDay = datetime.strptime(sArchiveEndDate, '%Y-%m-%d')
    except:
        wasdi.wasdiLog('End Date not valid, assuming today')

    oActualDate = oStartDay
    oTimeDelta = timedelta(days=1)
    oTimeDelta2 = oEndDay - oStartDay

    iDays = oTimeDelta2.days

    iStep = 100.0 / float(iDays + 1)
    iProgress = 0.0

    aoPayload = {}

    if sOrbits is None:
        sOrbits = ""

    if sOrbits == "":
        wasdi.wasdiLog("Orbits not found in input. Start Autosearch")
        sOrbits = findInvolvedOrbits(sArchiveStartDate, sArchiveEndDate, oBbox, sPROVIDER, iMaxSearchDays, sPlatformFilter)

    wasdi.wasdiLog("Saving the first version of the payload")
    aoPayload["inputs"] = wasdi.getParametersDict()
    aoPayload["orbits"] = sOrbits
    aoPayload["RUNS"] = []
    wasdi.setPayload(aoPayload)

    asWorkspaceFiles = wasdi.getProductsByActiveWorkspace()

    if sCopDemMap is None:
        sCopDemMap = ''

    if sCopDemMap == '':
        sCopDemMap = "__"

    if sLULCMap is None:
        sLULCMap = ''

    if sLULCMap == '':
        sLULCMap = "__"

    asLatLonSteps = sGridStep.split(",")
    dLatStep = float(asLatLonSteps[0])
    dLonStep = float(asLatLonSteps[1])

    oNewBbox = {}
    oNewBbox["southWest"] = {}
    oNewBbox["northEast"] = {}
    oNewBbox["southWest"]["lat"] = oBbox["southWest"]["lat"]
    oNewBbox["northEast"]["lat"] = oNewBbox["southWest"]["lat"]

    while oNewBbox["northEast"]["lat"] < oBbox["northEast"]["lat"]:
        oNewBbox["northEast"]["lat"] += dLatStep

    oNewBbox["southWest"]["lng"] = oBbox["southWest"]["lng"]
    oNewBbox["northEast"]["lng"] = oNewBbox["southWest"]["lng"]

    while oNewBbox["northEast"]["lng"] < oBbox["northEast"]["lng"]:
        oNewBbox["northEast"]["lng"] += dLonStep

    if sCopDemMap not in asWorkspaceFiles:
        wasdi.wasdiLog("Cop Dem Map Not found in input: try to generate it")

        aoCopDemFiles = wasdi.searchEOImages("StaticFiles", "2020-01-01", "2020-01-02", oBoundingBox=oNewBbox, sProductType="CopDEM30m_wbm_global_mosaic")

        if len(aoCopDemFiles)>0:
            wasdi.importProduct(aoCopDemFiles[0])
            sCopDemMap = aoCopDemFiles[0]["fileName"]
            wasdi.wasdiLog("Imported " + sCopDemMap + " as CopDem File")

    aoPayload["CopDemMap"] = sCopDemMap
    wasdi.setPayload(aoPayload)

    if not bUseHasard:
        if sLULCMap not in asWorkspaceFiles:
            wasdi.wasdiLog("LULC Map Not found in input: try to generate it")
            oLULCParams = {}
            sLULCMap = "LULC_Autowade_input_map.tif"
            oLULCParams["OUTPUT"] = sLULCMap
            oLULCParams["BBOX"] = oNewBbox
            oLULCParams["EXTRACT_MAP"] = False

            sLULCProcId = wasdi.executeProcessor("world_cover_extractor", oLULCParams)
            wasdi.waitProcess(sLULCProcId)

            wasdi.wasdiLog("Imported " + sLULCMap + " as World Cover File")

        aoPayload["LULC_IMAGE"] = sLULCMap

    while oActualDate <= oEndDay:
        sDate = oActualDate.strftime("%Y-%m-%d")

        wasdi.wasdiLog('-------------STARTING CHAIN FOR DATE ' + sDate)

        aoChainParams = {}
        aoChainParams["ENDDATE"] = sDate
        aoChainParams["DELETE"] = bDelete
        aoChainParams["SIMULATE"] = bSimulate
        aoChainParams["FORCE_RE_RUN"] = bForceReRun
        aoChainParams["BBOX"] = oBbox
        aoChainParams["ORBITS"] = sOrbits
        aoChainParams["GRIDSTEP"] = sGridStep
        aoChainParams["LASTDAYS"] = iLastDays
        aoChainParams["PREPROCWORKFLOW"] = sWorkflow
        aoChainParams["MOSAICBASENAME"] = sMosaicFileNameBase
        aoChainParams["SUFFIX"] = sSuffix
        aoChainParams["DETECTION_PROCESSOR"] = sFloodDetectionProcessor
        aoChainParams["LULC_IMAGE"] = sLULCMap

        # level depth of HSBA, default value is 8
        aoChainParams["depth"] = wasdi.getParameter("depth", 8)
        # Ashman D (AD) coefficient, default value is 2.6
        aoChainParams["ad"] = wasdi.getParameter("ad", 2.6)
        # initial threshold for singe image based water detection, default value is -17 dB
        aoChainParams["inithre"] = wasdi.getParameter("inithre", -17)
        # threshold of hierarchical image patch size of HSBA, default value is 5000 pixels
        aoChainParams["sizethre"] = wasdi.getParameter("sizethre", 5000)
        # size of blob removement, default value is 20 pixels
        # we may have a blob size and in case we use it
        # or if not maybe a noise_reduction, that is converted to a blob size
        # If none is available we assume the 20 default
        iBlobSize = wasdi.getParameter("blobsize", None)

        if iBlobSize == "":
            iBlobSize = None

        if iBlobSize is None:
            iNoiseReduction = wasdi.getParameter("NOISE_REDUCTION", None)

            if iNoiseReduction == "":
                iNoiseReduction = None

            if iNoiseReduction is None:
                iBlobSize = 85
            else:
                iBlobSize = int(20.0 + float(iNoiseReduction)*1.3)

        aoChainParams["blobsize"] = iBlobSize
        # bin size of histogram, default value is 0.3 dB
        aoChainParams["bin_step"] = wasdi.getParameter("bin_step", 0.3)
        # name of previous flood probability
        aoChainParams["pre_prob"] = wasdi.getParameter("pre_prob", None)
        # SAR shadow mask
        aoChainParams["sar_shadow"] = wasdi.getParameter("sar_shadow", None)
        # SAR low backscattering mask
        aoChainParams["low_bsc"] = wasdi.getParameter("low_bsc", None)
        # SAR non-sensitivity mask
        aoChainParams["nonsi_mask"] = wasdi.getParameter("nonsi_mask", None)
        # HAND mask
        aoChainParams["hand_mask"] = wasdi.getParameter("hand_mask", None)
        # water mask from CopDEM
        aoChainParams["copdem_wm"] = sCopDemMap
        # Apply map conversion or not?
        aoChainParams['APPLYMAPCONVERSION'] = bApplyMapConversion
        # Min Percentage of an image to intersect the shapefile or the area of interest
        aoChainParams["MINCOVERAGE"] = iMinimumPercentage
        # Flag to know if we need to move the result to sftp server or not
        aoChainParams["MOVETOSFTP"] = bMoveToSftp

        aoChainParams['NODATAVALUE'] = iNODATAVALUE
        aoChainParams['INPUTIGNOREVALUE'] = iINPUTIGNOREVALUE
        aoChainParams['MOSAICNODATAVALUE'] = iMOSAICNODATAVALUE
        aoChainParams['MOSAICINPUTIGNOREVALUE'] = iMOSAICINPUTIGNOREVALUE
        aoChainParams['PROVIDER'] = sPROVIDER
        aoChainParams['MAXREFERENCETILEAGE'] = iMaxReferenceTileAge
        aoChainParams['SHAPEFILEMASK'] = sShapeMaskFile
        aoChainParams['APPLY_WATER_MAP'] = bApplyWaterMap
        aoChainParams['PERMANENT_WATER_MAP_NAME'] = sPermanentWaterMapName

        aoChainParams["SLOPE_IMAGE"] = wasdi.getParameter("SLOPE_IMAGE", None)
        aoChainParams["DEM_IMAGE"] = wasdi.getParameter("DEM_IMAGE", None)

        # AUTOWADE PARAMETERS:
        aoChainParams["MIN_CLUST_N"] = iMinCluster
        aoChainParams["FILTER_SIZE"] = iFilterSize

        # FILTER ON PLATFORM NAMES
        aoChainParams["PLATFORM_FILTER"] = sPlatformFilter

        #sOuputFile = sMosaicFileNameBase + "_" + sDate + "_flood.tif"
        sOuputFile = sMosaicFileNameBase + "_" + sDate + "_" + sSuffix

        if sOuputFile in asWorkspaceFiles and bForceReRun:
            wasdi.wasdiLog('Output already present, force re run true: delete it')
            wasdi.deleteProduct(sOuputFile)
            asWorkspaceFiles.remove(sOuputFile)

        if sOuputFile not in asWorkspaceFiles:

            sProcessId = wasdi.executeProcessor("autofloodchain2", aoChainParams)

            wasdi.wasdiLog('Chain started waiting for end')

            sChainStatus = wasdi.waitProcess(sProcessId)

            aoRunResult = {"DATE": sDate, "PROCID": sProcessId, "STATUS": sChainStatus}

            aoPayload["RUNS"].append(aoRunResult)
            wasdi.wasdiLog('Chain done for day ' + sDate + " with status " + sChainStatus)
        else:
            wasdi.wasdiLog('Output file ' + sOuputFile + ' already in the workspace, jump')

        iProgress = iProgress + iStep
        wasdi.updateProgressPerc(int(iProgress))

        oActualDate += oTimeDelta

    if bDelete:
        if bApplyWaterMap and bApplyMapConversion:
            deletePermanentWaterTiles(oNewBbox["northEast"]["lat"], oNewBbox["southWest"]["lng"], oNewBbox["southWest"]["lat"], oNewBbox["northEast"]["lng"], dLatStep, dLonStep, sMosaicFileNameBase)

    wasdi.updateStatus("DONE", 100)
    wasdi.setPayload(aoPayload)

if __name__ == '__main__':
    wasdi.init("./config.json")
    run()