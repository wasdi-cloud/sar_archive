import wasdi
import numpy
from datetime import datetime
from datetime import timedelta
from osgeo import gdal, osr, ogr

def deleteFile(bDelete, sFileName):
    if bDelete:
        wasdi.deleteProduct(sFileName)

def hasMinIntersection(oBbox, oImage, fMinimumPercentage):

    try:

        sSouth = str(oBbox["southWest"]["lat"])
        sWest = str(oBbox["southWest"]["lng"])
        sNorth = str(oBbox["northEast"]["lat"])
        sEast = str(oBbox["northEast"]["lng"])

        # Create the WKT of the Target Area
        sTargetAreaWkt = "POLYGON (( " + sWest + " " + sSouth + "," + sWest + " " + sNorth + "," + sEast + " " + sNorth + "," + sEast + " " + sSouth + "," + sWest + " " + sSouth + "))"
        oTargetFootprint = ogr.CreateGeometryFromWkt(sTargetAreaWkt)

        # Get the image footprint
        sWktFootprint = oImage["footprint"]
        oImageFootprint = ogr.CreateGeometryFromWkt(sWktFootprint)

        # Intersect with target footprint
        oIntersection = oImageFootprint.Intersection(oTargetFootprint)
        # Get the area of the intersection
        fIntersectionArea = oIntersection.GetArea()
        # Get the area of the image
        fImageArea = oImageFootprint.GetArea()
        # Compute the percentage
        fCoverPercentage = fIntersectionArea / fImageArea

        if fCoverPercentage>fMinimumPercentage:
            return [True, fCoverPercentage]
        else:
            return [False, fCoverPercentage]

    except:
        wasdi.wasdiLog("Exception calculating intersection. Return Safe True")
        return True

def filterOnShapeFile(sShapeMaskFile, aoFoundImages, fMinPercentage):
    sShapeFullPath = wasdi.getPath(sShapeMaskFile)
    oShape = ogr.Open(sShapeFullPath)
    layer = oShape.GetLayer()
    feature = layer.GetFeature(0)
    vectorGeometry = feature.GetGeometryRef()
    #fShapeArea = vectorGeometry.GetArea()

    aoIncludedImages = []

    for oImage in aoFoundImages:
        try:
            # Get the image footprint
            sFootprint = oImage["footprint"]
            oGeometry = ogr.CreateGeometryFromWkt(sFootprint)
            # Get the area (surface) of the image footprint
            fImageArea = oGeometry.GetArea()

            # Intersect image footprint with shapefile
            oIntersection = oGeometry.Intersection(vectorGeometry)
            # Get the area of the intersection
            fIntersectionArea = oIntersection.GetArea()
            # Compute the percentage w.r.t. the area of the image footprint
            fCoverPercentage = fIntersectionArea / fImageArea

            #if oGeometry.Intersect(vectorGeometry):
            if fCoverPercentage > fMinPercentage:
                aoIncludedImages.append(oImage)
        except Exception as oE:
            wasdi.wasdiLog('Error intersecting image ' + str(oE))

    return aoIncludedImages

def generateDailySarMosaic(sMosaicFileName, aoFoundImages, sOrbit, aoPayload):
    # let see the images we are going to mosaic
    asMosaicInputs = []
    # Read the workflow for preprocessing
    sWorkflow = wasdi.getParameter('PREPROCWORKFLOW', 'LISTSinglePreproc2')
    # Read the data provider
    sPROVIDER = wasdi.getParameter('PROVIDER', 'AUTO')
    # Flag to know if delete is active or not
    bDelete = wasdi.getParameter('DELETE', False)
    # No Data Value
    iNODATAVALUE = wasdi.getParameter('NODATAVALUE', -9999)
    # Ignore Value
    iINPUTIGNOREVALUE = wasdi.getParameter('INPUTIGNOREVALUE', 0)

    # Log the number of images found for this orbit
    wasdi.wasdiLog('Found ' + str(len(aoFoundImages)) + ' images')
    aoPayload["ResultsPerOrbit"].append({"orbit": sOrbit, "images": len(aoFoundImages)})

    aoImagesToProcess = []
    asWorkspaceFiles = wasdi.getProductsByActiveWorkspace()

    for oImage in aoFoundImages:
        sFileName = oImage["fileName"]
        sPreprocessedFileName = sFileName.replace(".zip", "_preproc.tif")
        if sPreprocessedFileName not in asWorkspaceFiles:
            aoImagesToProcess.append(oImage)

    if len(aoImagesToProcess) > 0:
        # Import and pre-process all the images
        wasdi.importAndPreprocess(aoImagesToProcess, sWorkflow, '_preproc.tif', sPROVIDER)
        wasdi.wasdiLog('Workflows done. Clean original S1 files')
    else:
        wasdi.wasdiLog('All images are already imported and pre-processed')

    # Second cycle on images to find the ones to mosaic
    for iImages in range(0, len(aoFoundImages)):
        # Get the name
        sFileName = aoFoundImages[iImages]['fileName']
        # Create output name
        sPreprocessedName = sFileName.replace('.zip', '_preproc.tif')

        # Check if the preprocessed image is in place
        bExistsCheck = wasdi.fileExistsOnWasdi(sPreprocessedName)

        if bExistsCheck:
            aoPayload["InputImages"].append(sFileName)
            wasdi.wasdiLog('add ' + sPreprocessedName + ' to mosaic and delete original img [' + sFileName + ']')
            deleteFile(bDelete, sFileName)

            # Add the file for the mosaic
            asMosaicInputs.append(sPreprocessedName)
        else:
            aoPayload["FailedInputs"].append(sFileName)
            wasdi.wasdiLog('NOT ADDING File = ' + sPreprocessedName + ' file result not existing')

    iImagesToMosaic = len(asMosaicInputs)

    wasdi.wasdiLog('Running mosaic of ' + str(iImagesToMosaic) + ' images in file ' + sMosaicFileName)
    sStatus = wasdi.mosaic(asMosaicInputs, sMosaicFileName, int(iNODATAVALUE), int(iINPUTIGNOREVALUE))
    wasdi.wasdiLog('Mosaic done with status ' + sStatus)

    return asMosaicInputs

def filterImagesPerOrbit(aoAllImages, sOrbit):
    aoFilteredImages = []

    try:
        iRelativeOrbit = int(sOrbit)

        for oImage in aoAllImages:
            iActualImageOrbit = oImage["properties"]["relativeorbitnumber"]
            if int(iActualImageOrbit) == int(iRelativeOrbit):
                aoFilteredImages.append(oImage)

    except Exception as oE:
        wasdi.wasdiLog("filterImagesPerOrbit Error " + str(oE))

    return aoFilteredImages

def startFloodDetection(sReferenceImageName, sTile, asFloodProcessesToWait, asDemSubsetToDelete, asUncertantyToDelete, sUnderscoredSuffix):

    # Flag to force re-run
    bForceReRun = wasdi.getParameter('FORCE_RE_RUN', False)
    # Read Simulate Flag
    bSimulate = wasdi.getParameter('SIMULATE', False)

    # Log the reference
    wasdi.wasdiLog('Found reference image [' + sReferenceImageName + '] for tile [' + sTile + ']')

    # Generate names of reference flood map and output map
    #sReferenceFloodMap = sReferenceImageName.replace('.tif', '_flood.tif')
    sReferenceFloodMap = sReferenceImageName.replace('.tif', sUnderscoredSuffix)
    #sOutputFloodMap = sTile.replace('.tif', '_flood.tif')
    sOutputFloodMap = sTile.replace('.tif', sUnderscoredSuffix)

    # Check if the output is already available
    bOutputFloodMapExists = wasdi.fileExistsOnWasdi(sOutputFloodMap)

    # Should we force a delete?
    if bOutputFloodMapExists:
        if bForceReRun:
            wasdi.deleteProduct(sOutputFloodMap)
            bOutputFloodMapExists = False

    if not bOutputFloodMapExists:
        # SCHEDULE THE FLOOD PROCESSOR
        aoParams = {}

        # Name of reference image for change detection
        aoParams['ref_img'] = sReferenceImageName

        bLastFloodExists = wasdi.fileExistsOnWasdi(sReferenceFloodMap)
        if bLastFloodExists:
            wasdi.wasdiLog('Previous flood map found ' + sReferenceFloodMap)
            aoParams['pre_fld'] = sReferenceFloodMap

        # name of flood image for change detection
        aoParams["fld_img"] = sTile
        # Name of the output file Map
        aoParams["output_floodmap"] = sOutputFloodMap
        # level depth of HSBA, default value is 8
        aoParams["depth"] = wasdi.getParameter("depth", 8)
        # Ashman D (AD) coefficient, default value is 2.6
        aoParams["ad"] = wasdi.getParameter("ad", 2.6)
        # initial threshold for singe image based water detection, default value is -17 dB
        aoParams["inithre"] = wasdi.getParameter("inithre", -17)
        # threshold of hierarchical image patch size of HSBA, default value is 5000 pixels
        aoParams["sizethre"] = wasdi.getParameter("sizethre", 5000)
        # size of blob removement, default value is 20 pixels
        aoParams["blobsize"] = wasdi.getParameter("blobsize", 20)
        # bin size of histogram, default value is 0.3 dB
        aoParams["bin_step"] = wasdi.getParameter("bin_step", 0.3)
        # name of previous flood probability
        aoParams["pre_prob"] = wasdi.getParameter("pre_prob", None)
        # SAR shadow mask
        aoParams["sar_shadow"] = wasdi.getParameter("sar_shadow", None)
        # SAR low backscattering mask
        aoParams["low_bsc"] = wasdi.getParameter("low_bsc", None)
        # SAR non-sensitivity mask
        aoParams["nonsi_mask"] = wasdi.getParameter("nonsi_mask", None)
        # HAND mask
        aoParams["hand_mask"] = wasdi.getParameter("hand_mask", None)
        # CopDEM
        aoParams["copdem_wm"] = wasdi.getParameter("copdem_wm", None)

        if aoParams["pre_prob"] is None:
            aoParams["pre_prob"] = ""

        if aoParams["pre_prob"] == "":
            sPreProbMapName = "Uncertainty_" + sReferenceImageName
            aoParams["pre_prob"] = sPreProbMapName
            asUncertantyToDelete.append(sPreProbMapName)

        if not aoParams["copdem_wm"] is None:
            if aoParams["copdem_wm"] != '':
                wasdi.wasdiLog("Copdem available: warp on actual tile")
                sCopDem = warpMapOnTile(aoParams["copdem_wm"], sTile)
                wasdi.addFileToWASDI(sCopDem)
                asDemSubsetToDelete.append(sCopDem)
                aoParams["copdem_wm"] = sCopDem

        if not bSimulate:
            sFloodProcId = wasdi.asynchExecuteProcessor('hsba_pip', aoParams)
            wasdi.wasdiLog('Flood detection Processor scheduled with proc ID ' + sFloodProcId)
            asFloodProcessesToWait.append(sFloodProcId)
        else:
            wasdi.wasdiLog('Simulation mode on: fake scheduling')
    else:
        wasdi.wasdiLog('Flood tile already exists ' + sOutputFloodMap)

    return sOutputFloodMap


def startAUTOWADEFloodDetection(sReferenceImageName, sTile, asFloodProcessesToWait, asDemSubsetToDelete, asUncertantyToDelete, sUnderscoredSuffix):

    # Flag to force re-run
    bForceReRun = wasdi.getParameter('FORCE_RE_RUN', False)
    # Read Simulate Flag
    bSimulate = wasdi.getParameter('SIMULATE', False)

    # Log the reference
    wasdi.wasdiLog('Found reference image [' + sReferenceImageName + '] for tile [' + sTile + ']')

    # Generate names of reference flood map and output map
    sReferenceFloodMap = sReferenceImageName.replace('.tif', sUnderscoredSuffix)
    sOutputFloodMap = sTile.replace('.tif', sUnderscoredSuffix)

    # Check if the output is already available
    bOutputFloodMapExists = wasdi.fileExistsOnWasdi(sOutputFloodMap)

    # Should we force a delete?
    if bOutputFloodMapExists:
        if bForceReRun:
            wasdi.deleteProduct(sOutputFloodMap)
            bOutputFloodMapExists = False

    if not bOutputFloodMapExists:
        # SCHEDULE THE FLOOD PROCESSOR
        aoParams = {}

        # Name of reference image for change detection
        aoParams['PRE_IMAGE'] = sReferenceImageName
        # name of flood image for change detection
        aoParams["POST_IMAGE"] = sTile
        # Name of the output file Map
        aoParams["OUTPUT"] = sOutputFloodMap

        # level depth of HSBA, default value is 8
        aoParams["MIN_CLUST_N"] = wasdi.getParameter("MIN_CLUST_N", 6)
        # Ashman D (AD) coefficient, default value is 2.6
        aoParams["FILTER_SIZE"] = wasdi.getParameter("FILTER_SIZE", 3)
        aoParams["LULC_IMAGE"] = wasdi.getParameter("LULC_IMAGE", None)

        aoParams["SLOPE_IMAGE"] = wasdi.getParameter("SLOPE_IMAGE", None)
        aoParams["DEM_IMAGE"] = wasdi.getParameter("DEM_IMAGE", None)

        if not aoParams["LULC_IMAGE"] is None:
            if aoParams["LULC_IMAGE"] != '':
                wasdi.wasdiLog("LULC_IMAGE available: warp on actual tile")
                sCopDem = warpMapOnTile(aoParams["LULC_IMAGE"], sTile)
                wasdi.addFileToWASDI(sCopDem)
                asDemSubsetToDelete.append(sCopDem)
                aoParams["LULC_IMAGE"] = sCopDem

        if not bSimulate:
            sFloodProcId = wasdi.asynchExecuteProcessor('autowade_s1', aoParams)
            wasdi.wasdiLog('Flood detection Processor scheduled with proc ID ' + sFloodProcId)
            asFloodProcessesToWait.append(sFloodProcId)
        else:
            wasdi.wasdiLog('Simulation mode on: fake scheduling')
    else:
        wasdi.wasdiLog('Flood tile already exists ' + sOutputFloodMap)

    return sOutputFloodMap


def array2raster(oGeotransformation, oArray1, sRasterName='floodCountMap.tif', iEPSG=4326):
    """
    Function that converts a numpy array into a raster (.tif file) that can be easily visualized in QGIS.

    :param oGeotransformation: tuple with (originX, pixelWidth, 0, originY, 0, pixelHeight).
    From an exisiting 'data0' file, geotrans= data0.GetGeoTransform()
    The geotransform is used to convert from map to pixel coordinates and back using an affine transformation.
    The 3rd and 5th parameter are used (together with the 2nd and 4th) to define the rotation if your image doesn't have 'north up'.
    But most images are north up, and then both the 3rd and 5th parameter are zero.
    :param oArray1: 2D array of the first band to be converted into raster.
    :param sRasterName: the name of the raster file
    :param sEPSG: EPSG specification
    :return: Returns nothing, but stores the raster file directly to the specified directory.
    """
    try:
        iCols = oArray1.shape[1]
        iRows = oArray1.shape[0]

        bBigTiff = False
        if iCols * iRows * 8 > 4000000000:
            bBigTiff = True
            wasdi.wasdiLog('This looks like a big tiff')

        oDriver = gdal.GetDriverByName('GTiff')

        wasdi.wasdiLog(f'array2raster: create geotiff')

        iBands = 1

        aoOptions = ['COMPRESS=LZW']
        if bBigTiff:
            aoOptions.append('BIGTIFF=YES')

        oOutRaster = oDriver.Create(sRasterName, iCols, iRows, iBands, gdal.GDT_Byte, options=aoOptions)
        oOutRaster.SetGeoTransform(oGeotransformation)
        oOutBand1 = oOutRaster.GetRasterBand(1)
        oOutBand1.WriteArray(oArray1)

        oOutRasterSRS = osr.SpatialReference()
        oOutRasterSRS.ImportFromEPSG(iEPSG)
        oOutRaster.SetProjection(oOutRasterSRS.ExportToWkt())

        oOutBand1.FlushCache()

        wasdi.wasdiLog(f'array2raster: done :-)')
    except Exception as oE:
        wasdi.wasdiLog(f'array2raster: {type(oE)}: {oE}, aborting')


def warpMapOnTile(sMapToWarp, sTile, sOutputName=None):
    try:
        sMapPath = wasdi.getPath(sMapToWarp)
        sTilePath = wasdi.getPath(sTile)

        oTileDataset = gdal.Open(sTilePath)

        iRows = oTileDataset.RasterYSize
        iCols = oTileDataset.RasterXSize

        oTileGeoTransform = oTileDataset.GetGeoTransform()
        fPixelSizeX = oTileGeoTransform[1]
        fPixelSizeY = -oTileGeoTransform[5]

        dUpperLeftX = oTileGeoTransform[0]
        dUpperLeftY = oTileGeoTransform[3]

        dLowerRightX = dUpperLeftX + (iRows * fPixelSizeX)
        dLowerRightY = dUpperLeftY - (iCols * fPixelSizeY)

        oTileDataset = None

        if sOutputName is None:
            sOutputWarpedMap = sMapToWarp.replace(".tif", "_" + sTile)
        else:
            sOutputWarpedMap = sOutputName

        if not wasdi.fileExistsOnWasdi(sOutputWarpedMap):
            sWarpOptions = "-tr " + str(fPixelSizeX) + " " + str(fPixelSizeY)
            sWarpOptions = sWarpOptions + " "
            sWarpOptions = sWarpOptions + "-te " + str(dUpperLeftX) + " " + str(dLowerRightY) + " " + str(
                dLowerRightX) + " " + str(dUpperLeftY)
            sWarpOptions = sWarpOptions + "-of GTiff -co COMPRESS=DEFLATE "

            sOutputWarpedMapFullPath = wasdi.getPath(sOutputWarpedMap)
            gdal.Warp(sOutputWarpedMapFullPath, sMapPath, options=sWarpOptions)

        return sOutputWarpedMap

    except Exception as oE:
        wasdi.wasdiLog(f'warpWaterMapOnTile: {type(oE)}: {oE}, aborting')


def convertMapToFourOutputs(oBandAsArray, iNoData):

    iMax = numpy.max(oBandAsArray)

    if iMax>1 and iMax != iNoData:
        wasdi.wasdiLog("The array has values > 1, we suppose it is has been already converted")
        return oBandAsArray

    # Create the four state: nodata, noflood, permanent, flooded
    oBandAsArray = numpy.where(oBandAsArray == 1, 3, oBandAsArray)
    oBandAsArray = numpy.where(oBandAsArray == 0, 1, oBandAsArray)
    oBandAsArray = numpy.where(oBandAsArray == iNoData, 0, oBandAsArray)
    return oBandAsArray


def mosaicFloodMap(sFinalMapFile, asFloodMapsToMosaic, iMOSAICNODATAVALUE, iMOSAICINPUTIGNOREVALUE, bApplyMapConversion,
                   bApplyWaterMap, sPermanentWaterMapName, sUnderscoredSuffix, oBbox):
    try:
        sBaseName = wasdi.getParameter("MOSAICBASENAME", 'mosaic')

        # Divide the available maps per tile
        aoTilesPerOrbit = {}

        asFilesToDelete = []

        # For each map to mosaic
        for sFloodedTile in asFloodMapsToMosaic:
            # Slit the name parts
            asNameParts = sFloodedTile.split("_")
            if asNameParts is None:
                continue
            if len(asNameParts) < 6:
                continue
            if asNameParts[0] != sBaseName:
                continue

            # Get back only the x_y tile
            sTile = asNameParts[3] + "_" + asNameParts[4]

            # Did we already find this tile?
            if sTile not in aoTilesPerOrbit:
                aoTilesPerOrbit[sTile] = []

            # Append the image to the tile
            aoTilesPerOrbit[sTile].append(sFloodedTile)

        asFloodMapsToMosaic2 = []
        # Now we need to manually merge images if more than one is available for a single tile
        for sTile in aoTilesPerOrbit.keys():

            # How many images do we have for this tile?
            if len(aoTilesPerOrbit[sTile]) > 1:
                # More than one: so we will need to merge: create the final array
                aiCompositeData = []
                iAddedFiles = 0
                # Let name this as BASE_ALL_DATE_TILEX_TILEY + sUnderscoredSuffix
                asNameParts = sFloodedTile.split("_")
                sCumulatedTile = asNameParts[0] + "_ALL_" + asNameParts[2] + "_" + sTile + sUnderscoredSuffix
                oGeoTransform = None

                sFinalTileName = ""

                # For all the images of this tile
                for sFloodedTile in aoTilesPerOrbit[sTile]:
                    # Open the file
                    sFilePath = wasdi.getPath(sFloodedTile)
                    oInputDataset = gdal.Open(sFilePath)

                    # Read the band
                    iRows = oInputDataset.RasterYSize
                    iCols = oInputDataset.RasterXSize
                    oBand = oInputDataset.GetRasterBand(1)
                    oBandAsArray = oBand.ReadAsArray(0, 0, iCols, iRows)

                    if bApplyMapConversion:
                        # Create the four state: nodata, noflood, permanent, flooded
                        oBandAsArray = convertMapToFourOutputs(oBandAsArray, iMOSAICNODATAVALUE)

                        if iAddedFiles == 0:
                            # This is the first file: we can initialize the output array
                            aiCompositeData = numpy.array(oBandAsArray)
                            # And the output geotransformation
                            oGeoTransform = oInputDataset.GetGeoTransform()
                        else:
                            # Other image(s): the max is ok for the 3-state output
                            aiCompositeData = numpy.maximum(aiCompositeData, oBandAsArray)
                    else:
                        # In this case we stay with the 2-state map
                        if iAddedFiles == 0:
                            # This is the first file: we can initialize the output array
                            aiCompositeData = numpy.array(oBandAsArray)
                            # And the output geotransformation
                            oGeoTransform = oInputDataset.GetGeoTransform()
                        else:
                            # Other image(s): we need to keep the flooded pixels
                            aiCompositeData = numpy.where(oBandAsArray == 1, 1, aiCompositeData)
                            aiCompositeData = numpy.where((aiCompositeData == iMOSAICNODATAVALUE) & (oBandAsArray == 0), 0, aiCompositeData)

                    # Free the gdal resource
                    oInputDataset = None

                    # Keep count of the files added
                    iAddedFiles = iAddedFiles + 1
                    wasdi.wasdiLog("Flood tile " + sFloodedTile + " added to single day flood tile " + sCumulatedTile)

                # Ok this tile is done: save the new image of the tile
                array2raster(oGeoTransform, aiCompositeData, wasdi.getPath(sCumulatedTile))
                # Add to WASDI
                wasdi.addFileToWASDI(sCumulatedTile)
                # This will be the one added to the mosaic
                asFloodMapsToMosaic2.append(sCumulatedTile)
                # And remember to delete this temporary file (the ALL for this tile in this day)
                asFilesToDelete.append(sCumulatedTile)

                sFinalTileName = sCumulatedTile
            else:
                # There is only one image in this tile: so we just need to check if we must convert it
                sFloodedTile = aoTilesPerOrbit[sTile][0]

                if bApplyMapConversion:
                    # Yes, we need to convert
                    asNameParts = sFloodedTile.split("_")
                    # Prepare the ALL file name: the original tile must keep the binary values to let the chain work
                    #sCumulatedTile = asNameParts[0] + "_ALL_" + asNameParts[2] + "_" + sTile + "_flood.tif"
                    sCumulatedTile = asNameParts[0] + "_ALL_" + asNameParts[2] + "_" + sTile + sUnderscoredSuffix

                    sFinalTileName = sCumulatedTile

                    # Open the original tile
                    sFilePath = wasdi.getPath(sFloodedTile)
                    oInputDataset = gdal.Open(sFilePath)
                    iRows = oInputDataset.RasterYSize
                    iCols = oInputDataset.RasterXSize

                    # Read the band
                    oBand = oInputDataset.GetRasterBand(1)
                    oGeoTransform = oInputDataset.GetGeoTransform()
                    oBandAsArray = oBand.ReadAsArray(0, 0, iCols, iRows)

                    # Create the four state: nodata, noflood, permanent, flooded
                    oBandAsArray = convertMapToFourOutputs(oBandAsArray, iMOSAICNODATAVALUE)

                    # Save to a new raster, the ALL tile for this day
                    array2raster(oGeoTransform, oBandAsArray, wasdi.getPath(sCumulatedTile))

                    # Close the orginal file
                    oInputDataset = None

                    # Add the converted file to WASDI
                    wasdi.addFileToWASDI(sCumulatedTile)

                    # And remember to delete this temporary file (the ALL for this tile in this day)
                    asFilesToDelete.append(sCumulatedTile)

                    wasdi.wasdiLog("Flood tile " + sFloodedTile + " converted")
                else:
                    # No map conversion: we just need to mosaic the tile
                    sFinalTileName = sFloodedTile

                # We add the file to the list of files to mosaic
                asFloodMapsToMosaic2.append(sFinalTileName)

            if bApplyWaterMap and bApplyMapConversion:
                applyPermanetWaterMapToTile(sFinalTileName)

        # AOI GRID
        dLatN = oBbox["northEast"]["lat"]
        dLonW = oBbox["southWest"]["lng"]
        dLatS = oBbox["southWest"]["lat"]
        dLonE = oBbox["northEast"]["lng"]
        afBbox = [dLonW, dLatS, dLonE, dLatN]

        sFinalMapFileTmp = sFinalMapFile.replace(".tif", "_tmp.tif")

        if bApplyMapConversion:
            wasdi.wasdiLog("Saving daily 3-state flood map")
            # In this case 0 indicates nodata
            wasdi.mosaic(asFloodMapsToMosaic2, sFinalMapFileTmp, 0, 0)

            wasdi.wasdiLog("Cropping " + sFinalMapFileTmp + " to the bbox: " + str(afBbox))

            # to make sure to crop it to the bbox
            aoWarpOptions = gdal.WarpOptions(outputBounds=afBbox, creationOptions=["COMPRESS=LZW", "BIGTIFF=YES"], srcNodata=0, dstNodata=0)
            # arguments: name of the warped dataset (output), dataset to be warped (input), options
            gdal.Warp(wasdi.getPath(sFinalMapFile), wasdi.getPath(sFinalMapFileTmp), options=aoWarpOptions)
            wasdi.wasdiLog("Generated daily 3-state flood map: " + sFinalMapFile)

        else:
            wasdi.wasdiLog("Saving daily 2-state flood map")
            # In this case no data is the one of the mosaic, i.e. 255
            wasdi.mosaic(asFloodMapsToMosaic2, sFinalMapFileTmp, iMOSAICNODATAVALUE, iMOSAICNODATAVALUE)

            wasdi.wasdiLog("Cropping " + sFinalMapFileTmp + " to the bbox: " + str(afBbox))

            # to make sure to crop it to the bbox
            aoWarpOptions = gdal.WarpOptions(outputBounds=afBbox, creationOptions=["COMPRESS=LZW", "BIGTIFF=YES"], srcNodata=iMOSAICNODATAVALUE, dstNodata=iMOSAICNODATAVALUE)
            # arguments: name of the warped dataset (output), dataset to be warped (input), options
            gdal.Warp(wasdi.getPath(sFinalMapFile), wasdi.getPath(sFinalMapFileTmp), options=aoWarpOptions)
            wasdi.wasdiLog("Generated daily 2-state flood map: " + sFinalMapFile)

        # Add the map to WASDI
        wasdi.addFileToWASDI(sFinalMapFile)
        # delete the tmp file
        wasdi.deleteProduct(sFinalMapFileTmp)

        wasdi.wasdiLog("Flood mosaic " + sFinalMapFile + " done")

        # Here we clean the tile ALL images if present
        for sFile in asFilesToDelete:
            deleteFile(wasdi.getParameter("DELETE"), sFile)

    except Exception as oE:
        wasdi.wasdiLog("mosaicFloodMap Error " + str(oE))

    return

def getPermanentWaterMapNameFromTile(sInputTile):
    try:
        asNameParts = sInputTile.split("_")
        # Get the lat and lon index
        iTileLatIndex = int(asNameParts[3])
        iTileLonIndex = int(asNameParts[4])

        sPermanentWaterTileName = asNameParts[0]+"PERMWATER" + "_PW__" + str(iTileLatIndex) + "_" + str(iTileLonIndex) + ".tif"

        return sPermanentWaterTileName

    except Exception as oE:
        wasdi.wasdiLog("generatePermanentWaterTile Error " + str(oE))

    return ""
def generatePermanentWaterTile(sInputTile):
    try:
        asNameParts = sInputTile.split("_")
        # Get the lat and lon index
        iTileLatIndex = int(asNameParts[3])
        iTileLonIndex = int(asNameParts[4])

        # Total Bounding Box
        oBbox = wasdi.getParameter("BBOX")
        # Grid Lat,Lon Step
        sGridStep = wasdi.getParameter('GRIDSTEP', '1,1')
        # Get Lat Lon grid steps
        asLatLonSteps = sGridStep.split(",")

        # AOI GRID
        dLatN = oBbox["northEast"]["lat"]
        dLonW = oBbox["southWest"]["lng"]
        dLatS = oBbox["southWest"]["lat"]
        dLonE = oBbox["northEast"]["lng"]

        # GRID STEPS
        dLatStep = float(asLatLonSteps[0])
        dLonStep = float(asLatLonSteps[1])

        oTileBbox = {}
        oTileBbox["northEast"] = {}
        oTileBbox["southWest"] = {}

        oTileBbox["northEast"]["lat"] = dLatN - iTileLatIndex*dLatStep
        oTileBbox["southWest"]["lng"] = dLonW + iTileLonIndex * dLonStep
        oTileBbox["southWest"]["lat"] = dLatN - (iTileLatIndex+1) * dLatStep
        oTileBbox["northEast"]["lng"] = dLonW + (iTileLonIndex+1) * dLatStep

        aoWaterMapFiles = wasdi.searchEOImages("StaticFiles", "2020-01-01", "2020-01-02", oBoundingBox=oTileBbox, sProductType="ESACCI-Ocean-Land-Map-150m-P13Y-2000")

        if len(aoWaterMapFiles)>0:
            wasdi.importProduct(aoWaterMapFiles[0])
            sWaterMap = aoWaterMapFiles[0]["fileName"]
            wasdi.wasdiLog("Imported " + sWaterMap + " as Permanent Water File for tile " + str(iTileLatIndex) + "_" + str(iTileLonIndex))

            sPermanentWaterTileName = getPermanentWaterMapNameFromTile(sInputTile)
            warpMapOnTile(sWaterMap, sInputTile, sPermanentWaterTileName)

            wasdi.addFileToWASDI(sPermanentWaterTileName)

    except Exception as oE:
        wasdi.wasdiLog("generatePermanentWaterTile Error " + str(oE))

def applyPermanetWaterMapToTile(sTileName):

    sWaterMap = getPermanentWaterMapNameFromTile(sTileName)

    if sWaterMap == "":
        wasdi.wasdiLog("applyPermanetWaterMapToTile: error getting the permanet water map tile name")
        return

    asWorkspaceFiles = wasdi.getProductsByActiveWorkspace()

    if sWaterMap not in asWorkspaceFiles:
        generatePermanentWaterTile(sTileName)

    asWorkspaceFiles = wasdi.getProductsByActiveWorkspace()
    if sWaterMap in asWorkspaceFiles:
        sWaterMapTileFilePath = wasdi.getPath(sWaterMap)
        oWaterMapDataset = gdal.Open(sWaterMapTileFilePath)

        sTileFilePath = wasdi.getPath(sTileName)
        oTileDataset = gdal.Open(sTileFilePath, gdal.GA_Update)

        # Read the band
        iRows = oTileDataset.RasterYSize
        iCols = oTileDataset.RasterXSize
        oTileBand = oTileDataset.GetRasterBand(1)
        oTileBandAsArray = oTileBand.ReadAsArray(0, 0, iCols, iRows)

        oWaterMapBand = oWaterMapDataset.GetRasterBand(1)
        oWaterMapBandAsArray = oWaterMapBand.ReadAsArray(0, 0, iCols, iRows)

        oTileBandAsArray = numpy.where(oWaterMapBandAsArray == 2, 2, oTileBandAsArray)
        oTileBandAsArray = numpy.where(oWaterMapBandAsArray == 0, 2, oTileBandAsArray)

        oTileBand.WriteArray(oTileBandAsArray)
        oTileBand.FlushCache()

        oWaterMapDataset = None
        oTileDataset = None

    else:
        wasdi.wasdiLog("Impossible to find the water map " + sWaterMap + " in the workspace")
        return

def run():

    wasdi.wasdiLog("Auto Flood Chain 3.4.4")

    aoPayload = {
        "input": wasdi.getParametersDict()
    }
    wasdi.setPayload(aoPayload)

    # Read Bounding Box
    oBbox = wasdi.getParameter('BBOX')
    # Read Orbits
    sOrbits = wasdi.getParameter('ORBITS')
    # Read Search range in days
    iLastDays = wasdi.getParameter('LASTDAYS', 0)
    # Read Base Mosaic Name
    sMosaicFileNameBase = wasdi.getParameter('MOSAICBASENAME', 'mosaic')
    # Read suffix to apply to the daily flood maps
    sSuffix = wasdi.getParameter('SUFFIX', 'baresoil-flood.tif')
    # Grid Lat,Lon Step
    sGridStep = wasdi.getParameter('GRIDSTEP', '1,1')
    # Read End Date
    sEndDate = wasdi.getParameter('ENDDATE')
    # Delete flag
    bDelete = wasdi.getParameter('DELETE', False)
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
    # Flag to know if we need to apply the three state map conversion or not
    bApplyMapConversion = wasdi.getParameter("APPLYMAPCONVERSION", False)
    # Flag to know if we need to move the result to sftp server or not
    bMoveToSftp = wasdi.getParameter("MOVETOSFTP", False)
    # Shape file mask
    sShapeMaskFile = wasdi.getParameter('SHAPEFILEMASK')
    # Min requested cover of the AOI for an image to be considered
    iMinimumPercentage = wasdi.getParameter('MINCOVERAGE', 10)
    # Flag to know if we need to apply a permanent water map
    bApplyWaterMap = wasdi.getParameter('APPLY_WATER_MAP', False)
    # Name of the Permanent Water Map Name. If none will take the copernicus one
    sPermanentWaterMapName = wasdi.getParameter('PERMANENT_WATER_MAP_NAME', "")

    # Flood Detection
    sFloodDetectionProcessor = wasdi.getParameter("DETECTION_PROCESSOR", "HASARD")

    # Filter on Platform names: S1A, S1B, S2C
    sPlatformFilter = wasdi.getParameter("PLATFORM_FILTER", "")

    if sFloodDetectionProcessor not in ["HASARD", "AUTOWADE"]:
        sFloodDetectionProcessor = "HASARD"

    bUseHasard = True

    if sFloodDetectionProcessor == "AUTOWADE":
        bUseHasard = False
        bApplyMapConversion = True
        wasdi.wasdiLog("Using Autowade as Detection Processor: forcing Apply Map Conversion = True")

    if not sSuffix.startswith("_"):
        sUnderscoredSuffix = '_' + sSuffix

    if bApplyWaterMap is True and bApplyMapConversion is False:
        wasdi.wasdiLog("A water map is required, but Apply Map conversion is false. Impossible to add permanent water when producing binary maps: permanent water will be ignored")

    # Split the involved orbits
    asOrbits = sOrbits.split(",")

    # Get Lat Lon grid steps
    asLatLonSteps = sGridStep.split(",")

    # Get Time Interval
    iLastDays = int(iLastDays)

    # Init dates
    iNow = datetime.today()

    # Cast dates from params
    try:
        iNow = datetime.strptime(sEndDate, '%Y-%m-%d')
    except:
        wasdi.wasdiLog('End date not set: it will be set to today')

    if iLastDays == -1:
        wasdi.wasdiLog('iLastDays = -1 => It will be set to yesterday (UTC)')
        oTimeDeltaBack = timedelta(days=1)
        iNow = iNow - oTimeDeltaBack
        iLastDays = 0

    # Find the start date
    iStart = iNow - timedelta(days=iLastDays)

    sStartDate = iStart.strftime("%Y-%m-%d")
    sEndDate = iNow.strftime("%Y-%m-%d")

    wasdi.wasdiLog('Start Date = ' + sStartDate + ' End Date = ' + sEndDate)

    aoPayload["input"] = wasdi.getParametersDict()

    #Orbits
    aoPayload["Orbits"] = sOrbits
    #Input Images
    aoPayload["InputImages"] = []
    aoPayload["FailedInputs"] = []
    aoPayload["ResultsPerOrbit"] = []
    aoPayload["OutputAlreadyDone"] = False
    aoPayload["MosaicAlreadyDone"] = False

    sFinalMapFile = sMosaicFileNameBase + '_' + sEndDate + sUnderscoredSuffix

    bFinalMapAlreadyAvailable = wasdi.fileExistsOnWasdi(sFinalMapFile)

    if bForceReRun:
        if bFinalMapAlreadyAvailable:
            wasdi.wasdiLog("Force Re Run True: delete existing output map")
            wasdi.deleteProduct(sFinalMapFile)
            bFinalMapAlreadyAvailable = False

    if bFinalMapAlreadyAvailable:
        aoPayload["OutputAlreadyDone"] = True
        aoPayload["MosaicAlreadyDone"] = True
        aoPayload["OutputFile"] = sFinalMapFile
        wasdi.setPayload(aoPayload)
        wasdi.updateStatus("DONE", 100)
        wasdi.wasdiLog('Flood Map ' + sFinalMapFile + ' already exists: procedure Done')
        return

    asFloodProcessesToWait = []

    iPerc = 0

    fPercPerOrbit = 90.0

    if len(asOrbits) > 0:
        fPercPerOrbit = 90.0 / float(len(asOrbits))

    asTotalMosaicInputs = []
    asFloodMapsToMosaic = []
    asDemSubsetToDelete = []
    asUncertantyToDelete = []

    # Search all the images, for all the orbits
    #aoAllOrbitImages = wasdi.searchEOImages('S1', sStartDate, sEndDate, sProductType='GRD', iOrbitNumber=None, oBoundingBox=oBbox, sProvider=sPROVIDER)

    # if needed, filter based on platform name(s) for the search of all the images, for all the orbits
    if sPlatformFilter == "":
        wasdi.wasdiLog("No filter applied to the search of S1 images")
        aoAllOrbitImages = wasdi.searchEOImages('S1', sStartDate, sEndDate, sProductType='GRD', iOrbitNumber=None, oBoundingBox=oBbox, sProvider=sPROVIDER)
    else:
        wasdi.wasdiLog("Only searching for S1 images of platforms: " + str(sPlatformFilter))
        #we need to remove whitespaces!!!!!!!!!!!!!
        asPlatformFilter = [x.strip() for x in sPlatformFilter.split(',')]
        aoAllOrbitImages = []
        aoAllOrbitImagesAllPlatforms = wasdi.searchEOImages('S1', sStartDate, sEndDate, sProductType='GRD', iOrbitNumber=None, oBoundingBox=oBbox, sProvider=sPROVIDER)

        for oAllOrbitImageAllPlatforms in aoAllOrbitImagesAllPlatforms:
            if oAllOrbitImageAllPlatforms["fileName"].startswith(tuple(asPlatformFilter)):
                aoAllOrbitImages.append(oAllOrbitImageAllPlatforms)



    # If we do not have images, that is a problem for this day
    if len(aoAllOrbitImages) == 0:
        wasdi.wasdiLog('No images available in range ' + sStartDate + ' - ' + sEndDate)
        wasdi.setPayload(aoPayload)
        wasdi.updateStatus("DONE", 100)
        return

    fMinimumPercentage = 0.0
    try:
        if iMinimumPercentage is int:
            if iMinimumPercentage<0:
                iMinimumPercentage = 0
            if iMinimumPercentage>100:
                iMinimumPercentage = 100
            fMinimumPercentage = iMinimumPercentage/100
    except:
        wasdi.wasdiLog("Exception converting iMinimumPercentage to float")

    # Check if we have to filter input images based on a shape file
    bFilteredByShape = False
    try:
        if sShapeMaskFile is not None:
            if sShapeMaskFile != "":
                aoAllOrbitImages = filterOnShapeFile(sShapeMaskFile, aoAllOrbitImages, fMinimumPercentage)
                bFilteredByShape = True
    except:
        wasdi.wasdiLog("Exception accessing the shapefile. Geographic filter not active: the processor will take all the images in the workspace")

    # Check if we have to filter input images based on the percentual coverage of the Area of Interest
    if not bFilteredByShape and fMinimumPercentage > 0:
        try:
            aoFilteredAllOrbitImge = []
            for oImage in aoAllOrbitImages:
                [bMinIntersection, fCoverage] = hasMinIntersection(oBbox, oImage, fMinimumPercentage)
                if bMinIntersection:
                    aoFilteredAllOrbitImge.append(oImage)
                else:
                    iCoverage = int(fCoverage*100)
                    wasdi.wasdiLog("Image " + oImage["fileName"] + " does not meet min coverage (" + str(iCoverage) + "):" + " excluded")

            aoAllOrbitImages = aoFilteredAllOrbitImge
        except:
            wasdi.wasdiLog("Exception filtering with respect to min coverage")

    # For each orbit
    for sOrbit in asOrbits:
        # Generate the Mosaic File Name
        sMosaicFileName = sMosaicFileNameBase + '_' + sOrbit + '_' + sEndDate + '.tif'
        # Does it already exists?
        bMosaicExists = wasdi.fileExistsOnWasdi(sMosaicFileName)

        # If exists but we force re-run we have to delete it
        if bMosaicExists and bForceReRun:
            wasdi.deleteProduct(sMosaicFileName)
            bMosaicExists = False

        # Final Check if exists at the end or not
        if not bMosaicExists:
            # No we have to make it
            # Get the list of available images for this orbit
            aoFoundImages = filterImagesPerOrbit(aoAllOrbitImages, sOrbit)

            # Do we have any image?
            if len(aoFoundImages) == 0:
                continue

            wasdi.wasdiLog('Found images for orbit = ' + sOrbit)
            # Generate the mosaic
            asMosaicInputs = generateDailySarMosaic(sMosaicFileName, aoFoundImages, sOrbit, aoPayload)
            # Keep the list of images that will have to be deleted
            asTotalMosaicInputs.extend(asMosaicInputs)
        else:
            aoPayload["MosaicAlreadyDone"] = True
            wasdi.wasdiLog('Mosaic already exists for orbit ' + sOrbit)

        # AOI GRID
        dLatN = oBbox["northEast"]["lat"]
        dLonW = oBbox["southWest"]["lng"]
        dLatS = oBbox["southWest"]["lat"]
        dLonE = oBbox["northEast"]["lng"]

        # GRID STEPS
        dLatStep = float(asLatLonSteps[0])
        dLonStep = float(asLatLonSteps[1])

        iLatStepCount = 0
        iLonStepCount = 0

        #MOSAIC BBOX
        sMosaicBbox = wasdi.getProductBBOX(sMosaicFileName)

        if sMosaicBbox is None:
            wasdi.wasdiLog('Mosaic bbox is null: leaving bbox equal to AOI')
            sMosaicBbox = str(oBbox["southWest"]["lat"]) + "," + str(oBbox["southWest"]["lng"]) + ",," + str(oBbox["northEast"]["lng"]) + "," + str(oBbox["northEast"]["lat"])

        asMosaicBBoxPoints = sMosaicBbox.split(",")
        dMosaicLatN = float(asMosaicBBoxPoints[4])
        dMosaicLonW = float(asMosaicBBoxPoints[1])
        dMosaicLatS = float(asMosaicBBoxPoints[0])
        dMosaicLonE = float(asMosaicBBoxPoints[3])

        # MULTI Subset: create the array of subset to generate
        asGeneratedTiles = []
        asLatNList = []
        asLonWList = []
        asLatSList = []
        asLonEList = []

        for dLon in numpy.arange(dLonW, dLonE, dLonStep):
            for dLat in numpy.arange(dLatS, dLatN, dLatStep):
                sLatNorth = str(dLat + dLatStep)
                sLatSouth = str(dLat)
                sLonWest = str(dLon)
                sLonEast = str(dLon + dLonStep)

                sTileFileName = sMosaicFileNameBase + '_' + sOrbit + '_' + sEndDate + '_' + str(iLatStepCount) + '_' + str(iLonStepCount) + '.tif'

                if (float(sLatSouth) > dMosaicLatN) or (float(sLatNorth) < dMosaicLatS) or (float(sLonEast) < dMosaicLonW) or (float(sLonWest) > dMosaicLonE):
                    #wasdi.wasdiLog('tile ' + tileFileName + ' out of the mosaic bbox')
                    pass
                else:
                    bExistsCheck = wasdi.fileExistsOnWasdi(sTileFileName)

                    if bForceReRun:
                        wasdi.deleteProduct(sTileFileName)
                        bExistsCheck = False

                    if not bExistsCheck:
                        #wasdi.wasdiLog('Add Subset file ' + tileFileName + '  bbox LATN = ' + sLatNorth + ' LONW = ' + sLonWest + ' LATS = ' + sLatSouth + ' LONE = ' + sLonEast)
                        asGeneratedTiles.append(sTileFileName)
                        asLatNList.append(sLatNorth)
                        asLatSList.append(sLatSouth)
                        asLonEList.append(sLonEast)
                        asLonWList.append(sLonWest)
                    else:
                        wasdi.wasdiLog('tile ' + sTileFileName + ' already in place')

                iLatStepCount = iLatStepCount + 1

            iLatStepCount = 0

            iLonStepCount = iLonStepCount + 1

        if len(asGeneratedTiles) > 0:
            asActualTiles = []
            asActualLatN = []
            asActualLonW = []
            asActualLatS = []
            asActualLonE = []

            iLimitCounter = 1

            for iActualTileCount in range(0, len(asGeneratedTiles)):
                asActualTiles.append(asGeneratedTiles[iActualTileCount])
                asActualLatN.append(asLatNList[iActualTileCount])
                asActualLonW.append(asLonWList[iActualTileCount])
                asActualLatS.append(asLatSList[iActualTileCount])
                asActualLonE.append(asLonEList[iActualTileCount])

                if iLimitCounter == 10:
                    wasdi.wasdiLog('Running multiSubset of tiles# ' + str(len(asActualTiles)))
                    sStatus = wasdi.multiSubset(sMosaicFileName, asActualTiles, asActualLatN, asActualLonW, asActualLatS, asActualLonE)
                    wasdi.wasdiLog('Subset status ' + sStatus)
                    iLimitCounter = 0
                    asActualTiles = []
                    asActualLatN = []
                    asActualLonW = []
                    asActualLatS = []
                    asActualLonE = []

                iLimitCounter = iLimitCounter + 1

            if len(asActualTiles) > 0:
                sStatus = wasdi.multiSubset(sMosaicFileName, asActualTiles, asActualLatN, asActualLonW, asActualLatS, asActualLonE)
                wasdi.wasdiLog('Subset status ' + sStatus)
        else:
            wasdi.wasdiLog('No tiles to generate')

        # Array Of Generated Tiles
        asGeneratedTiles = []
        iLatStepCount = 0
        iLonStepCount = 0

        wasdi.wasdiLog('Check available tiles...')

        for dLon in numpy.arange(dLonW, dLonE, dLonStep):
            for dLat in numpy.arange(dLatS,  dLatN, dLatStep):
                sLatNorth = str(dLat + dLatStep)
                sLatSouth = str(dLat)
                sLonWest = str(dLon)
                sLonEast = str(dLon + dLonStep)

                sTileFileName = sMosaicFileNameBase + '_' + sOrbit + '_' + sEndDate + '_' + str(iLatStepCount) + '_' + str(iLonStepCount) + '.tif'

                # Check if the tile is fully out of the Mosaic BBOX
                if (float(sLatSouth) > dMosaicLatN) or (float(sLatNorth) < dMosaicLatS) or ( float(sLonEast) < dMosaicLonW) or (float(sLonWest) > dMosaicLonE):
                    # NOT MY TILE
                    #wasdi.wasdiLog("jump tile " + tileFileName)
                    pass
                else:
                    bExistsCheck = wasdi.fileExistsOnWasdi(sTileFileName)
                    if bExistsCheck:
                        asGeneratedTiles.append(sTileFileName)

                iLatStepCount = iLatStepCount + 1

            iLatStepCount = 0

            iLonStepCount = iLonStepCount + 1

        wasdi.wasdiLog('Available ' + str(len(asGeneratedTiles)) + ' tiles. Searching for reference images')

        # Get the list of products of theWorkspace(updated)
        asWorkspaceProducts = wasdi.getProductsByActiveWorkspace()

        # for each tile
        for iTiles in range(0, len(asGeneratedTiles)):
            #wasdi.wasdiLog('Searching reference image for tile ' + asGeneratedTiles[iTiles])

            # Split the name
            sTileNameWithoutExt = asGeneratedTiles[iTiles].replace('.tif', '')
            asTileNameParts = sTileNameWithoutExt.split("_")

            if asTileNameParts[0] != sMosaicFileNameBase:
                continue

            # Take the date part and remove - char
            sTileTimeStamp = asTileNameParts[2].replace('-', '')
            # Convert as a timestamp
            iTileTimeStamp = int(sTileTimeStamp)

            try:
                oTileDate = datetime.strptime(asTileNameParts[2], '%Y-%m-%d')
            except:
                wasdi.wasdiLog('Exception computing date of tile')

            # Initialize reference image
            iReferenceImageTimeStamp = 0
            sReferenceImageName = None

            asTilesToDelete = []

            # Search all the workspace
            for iFileNames in range(0, len(asWorkspaceProducts)):
                # Clean and split the file name
                sWSFileName = asWorkspaceProducts[iFileNames]
                sWSFileName = sWSFileName.replace('.tif', '')
                asFileNameParts = sWSFileName.split("_")

                if len(asFileNameParts) >= 5:

                    #Check Country
                    if asFileNameParts[0] != asTileNameParts[0]:
                        continue
                    #Check orbit
                    if asFileNameParts[1] != asTileNameParts[1]:
                        continue
                    #Check tile grid
                    if asFileNameParts[3] != asTileNameParts[3]:
                        continue
                    if asFileNameParts[4] != asTileNameParts[4]:
                        continue
                    # If it has 6 parts is not a Tile but a flood Map
                    if len(asFileNameParts) == 6:
                        continue
                    # is it myself?
                    if sTileNameWithoutExt == sWSFileName:
                        continue

                    # Ok found a good candidate

                    # Take its timestamp
                    sCandidateTimeStamp = asFileNameParts[2].replace('-', '')
                    iCandidateTimeStamp = int(sCandidateTimeStamp)

                    try:
                        oCandidateDate = datetime.strptime(asFileNameParts[2], '%Y-%m-%d')
                    except:
                        wasdi.wasdiLog('Error converting date of the candidate')

                    oTileCandidateTimeDelta = oTileDate-oCandidateDate

                    # It must be younger than the tile and older than all other candidates but not too old
                    if (iCandidateTimeStamp < iTileTimeStamp) and (iCandidateTimeStamp > iReferenceImageTimeStamp) and (oTileCandidateTimeDelta.days<iMaxReferenceTileAge):
                        iReferenceImageTimeStamp = iCandidateTimeStamp
                        #wasdi.wasdiLog('Set as reference image ' + sWSFileName)

                        if sReferenceImageName is not None:
                            #wasdi.wasdiLog('updated reference image: set to delete the old one ' + sReferenceImageName)
                            asTilesToDelete.append(sReferenceImageName)

                        sReferenceImageName = sWSFileName + '.tif'
                    else:
                        if iCandidateTimeStamp < iTileTimeStamp:
                            #wasdi.wasdiLog('candidate older than another one or too old: set to delete ' + asWorkspaceProducts[iFileNames])
                            asTilesToDelete.append(asWorkspaceProducts[iFileNames])

            # There was a reference image?
            if sReferenceImageName is not None:

                if bUseHasard:
                    # Start the real flood detection
                    sOutputFileMap = startFloodDetection(sReferenceImageName, asGeneratedTiles[iTiles], asFloodProcessesToWait, asDemSubsetToDelete, asUncertantyToDelete, sUnderscoredSuffix)
                else:
                    sOutputFileMap = startAUTOWADEFloodDetection(sReferenceImageName, asGeneratedTiles[iTiles], asFloodProcessesToWait, asDemSubsetToDelete, asUncertantyToDelete, sUnderscoredSuffix)

                # add it the list of daily files
                asFloodMapsToMosaic.append(sOutputFileMap)
            else:
                wasdi.wasdiLog('Reference image not found')

            if len(asTilesToDelete) > 0:
                wasdi.wasdiLog('Found tiles to delete: ' + str(len(asTilesToDelete)))
                for iTilesToDelete in range(0, len(asTilesToDelete)):
                    deleteFile(bDelete, asTilesToDelete[iTilesToDelete])
                    #deleteFile(bDelete, asTilesToDelete[iTilesToDelete].replace('.tif', '_flood.tif'))
                    deleteFile(bDelete, asTilesToDelete[iTilesToDelete].replace('.tif', sUnderscoredSuffix))

        wasdi.wasdiLog('Delete mosaic of pre-processed images: ' + sMosaicFileName)
        deleteFile(bDelete, sMosaicFileName)

        iPerc = iPerc + int(fPercPerOrbit)

        if iPerc > 90:
            iPerc = 90

        wasdi.updateProgressPerc(iPerc)

    wasdi.updateProgressPerc(90)

    if len(asFloodProcessesToWait) > 0:
        wasdi.wasdiLog('Waiting for all flood maps to finish ' + str(len(asFloodProcessesToWait)))
        asResults = wasdi.waitProcesses(asFloodProcessesToWait)

    if len(asDemSubsetToDelete) > 0:
        if bDelete:
            for sDemToDelete in asDemSubsetToDelete:
                wasdi.deleteProduct(sDemToDelete)

    if len(asUncertantyToDelete) > 0:
        if bDelete:
            for sUncertantyFile in asUncertantyToDelete:
                wasdi.deleteProduct(sUncertantyFile)

    if len(asFloodMapsToMosaic) > 0:
        # Schedule the mosaic of the flood map
        wasdi.wasdiLog('Creating flood map mosaic for date ' + sEndDate)
        mosaicFloodMap(sFinalMapFile, asFloodMapsToMosaic, iMOSAICNODATAVALUE, iMOSAICINPUTIGNOREVALUE, bApplyMapConversion, bApplyWaterMap, sPermanentWaterMapName, sUnderscoredSuffix, oBbox)

    if bMoveToSftp:
        wasdi.wasdiLog('Moving ' + sFinalMapFile + ' to sftp server folder')
        wasdi.copyFileToSftp(sFinalMapFile, bAsynch=True)

    wasdi.wasdiLog('Delete preprocessed images')
    # Clear  Pre - Processed Images
    for iImages in range(0, len(asTotalMosaicInputs)):
        deleteFile(bDelete, asTotalMosaicInputs[iImages])

    aoPayload["OutputFile"] = sFinalMapFile
    wasdi.setPayload(aoPayload)
    wasdi.updateStatus("DONE", 100)
    wasdi.wasdiLog('Done')

if __name__ == '__main__':
    wasdi.init("./config.json")
    run()
