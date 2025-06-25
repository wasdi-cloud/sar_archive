### Automatic Flood Detection Chain 2.1

This app runs an automatic process to search available S1 GRD Images over an area of interest, pre-process the images, create a mosaic for each available orbit and then triggers the flood detection for each tile obtained from the mosaics. The final output flood map is again mosaiced in a single flood map.

The chain is designed to produce a single flood map, usually for a single day. This processor can be called by edrift_archive_generator to produce the daily map for a full period.

The chain calls the hsba_pip processor for the "real" flood detection of single tiles.

The output map legend is:
0 = No Data
1 = Not-Flooded
2 = Permanent Water
3 = Flooded Area

#### Parameters

{  
**"ENDDATE":"2019-07-08"**, *//Reference Day of the map to produce in the syntax: YYYY-MM-DD*  
**"BBOX":"29.0,92.0,10.0,100.0"**, *// Bounding Box of the Area of Interest in the form: LATN,LONW,LATS,LONE*  
**"ORBITS":"33,41,62,70,77,99,106,135,143,172,4"**, *// Comma separated list of the S1 Orbits to consider*  
**"MOSAICBASENAME":"FL"**, *// Base Name of the output file. DO NOT USE "\_" or " " or other special chars in the Base Name* 
**"SUFFIX":"baresoil-flood.tif"**, *// it is a short string used as suffix for the name of the output maps. It has to be WITHOUT underscores and it must end in ".tif"*
**"DELETE":"1"**, *// Leave 1 to delete intermediate files. 0 to keep all intermediate files*  
**"SIMULATE":"0"**, *// Leave 0 for real execution. 1 does not trigger the real flood detection*  
**"PREPROCWORKFLOW":"LISTSinglePreproc2"**, *// Name of the SNAP Workflow used to pre-process S1 Images*  
**"HSBA\_DEPTH\_IN":"-1"**, *// Parameter of the edriftlistflood_archive processor: HSBA starting tilling level. -1 for automatic detection*  
**"ASHMAN\_COEFF":"2.4"**, *// Parameter of the edriftlistflood_archive processor:  Bimodal distribution check (Ashman coefficient)*  
**"MIN\_PIXNB\_BIMODD":"10000"**, *// Parameter of the edriftlistflood_archive processor: Minimum number of pixels per image tiles to be processed*  
**"BLOBS\_SIZE":"150"**, *// Parameter of the edriftlistflood_archive processor: Size (in pixels) of small objects to remove*  
**"GRIDSTEP":"1,1"**, *// LEGACY: leave as it is*  
**"LASTDAYS":"0"**,  *// LEGACY: Number of days to consider before end date. 0, default, consider only the end date (usually do not change this. To make a period run this chain for different days changing END Date)*  
**"MOSAICXSTEP":"-1.0"**, *// LEGACY: output mosaic x step. Leave -1 for default native resolution*    
**"MOSAICYSTEP":"-1.0"**, *// LEGACY: output mosaic x step. Leave -1 for default native resolution*  
**"NODATAVALUE":"-9999"**, *// ADVANCED: value to use as no data for the intermediate mosaic image*    
**"INPUTIGNOREVALUE":"0"**,  *// ADVANCED: value to ignore from the pre-processed S1 image*  
**"FLOODNODATAVALUE":"255"**, *// ADVANCED: value to use as no data for the output flood map tiles*  
**"FLOODINPUTIGNOREVALUE":"-9999"**, *// ADVANCED: value to ignore in input for the flood tile computation*    
**"MOSAICNODATAVALUE":"255"**, *// ADVANCED: value to use as no data value for output flood mosaic map*  
**"MOSAICINPUTIGNOREVALUE":"255"**, *// ADVANCED: value to ignore in input for the flood mosaic computation*      
**"APPLYMAPCONVERSION":"1"** *// ADVANCED: value to ignore in input for the flood mosaic computation*  
}