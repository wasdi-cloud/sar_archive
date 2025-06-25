## SAR Archive Generator 3.3.7

The SAR Archive Generator is used to generate daily flood maps in the time range and over the area of interest specified by the user. 
It detects floods in open areas, using Sentinel-1 SAR images. Maps can be generated as binary or 3 state (flooded, not-flooded, permanent water).
The basic usage is to select a start date (ARCHIVE_START_DATE), an end date (ARCHIVE_END_DATE) and the area of interest (BBOX).

The Base Name (BASENAME) is a short code used as prefix for the name of the output maps. It must be WITHOUT spaces or "_" underscores. The system will remove spaces or underscores if it finds them.
The SUFFIX is a short string used as suffix for the name of the output maps. It has to be WITHOUT underscores and it must end in ".tif".

More advanced options are available.

### Images Filter / Orbits Detection
The user can have some control to specify the images that will be used to create the archive.

* Orbits: the ORBITS paramter is a string with a comma separated list of the S1 relative orbits to include. It the parameter is left blank, the processor will perfom an automatic search and selection of the orbits to use. In this case, the system will search all the orbits that interesects the given area, starting from end date to start date with the constraint that will search MAX_SEARCH_DAYS_FOR_ORBITS days before maximum. A manual specification of the orbits can avoid to process useless images and obtain a faster result when possible
* Min Coverage Area: (MINCOVERAGE) defines the minimum percentage of coverage that an image must have, intersecting the Area of Interest, to be included. By default is zero; this means that even an image that cover 1% of our area of interest is considered. Increase this value to avoid to process images that does not give a huge value added
* Shape File Mask: (SHAPEFILEMASK) if you have a shape file that better defines the area of interest, you can upload it in a workspace, run the archive on the selected workspace, and specify here the name of the shape file. Only the images that interesects the shape file in the bounding box will be included

### Outputs
The processor generates a daily flooded area map for each day where a pre and post S1 GRD Image is found in the area of interest, starting from start date up to end date.
The output maps will have names like:

    CODE_YYYY-MM-DD_flood.tif

If the flag APPLYMAPCONVERSION is set to *false* the output maps are binary:

    0 = Not Flooded
    1 = Flooded
    (255 = No Data)

If the flag APPLYMAPCONVERSION is set to *true* the output maps are:

    (0 = No Data)
    1 = Not Flooded
    2 = Permanent Water
    3 = Flooded

### Permanent Water
If the APPLYMAPCONVERSION is set to true the system can overlay permanent water layer to the outputs. To do it select the "Apply Permanent Water" (APPLY_WATER_MAP) flag in the advanced Tab. The user can have a personal Permanent Water Map: in this case just create a Workspace, upload the Permanent Water Map, Open the application, select the workspace and select the upladed file in the Advanced Tab -> Permanent Water option.
If the APPLYMAPCONVERSION flag is true, the APPLY_WATER_MAP is true but no map is selected, the application will automatically extract the Permanent Water Map from the ESACCI-Ocean-Land-Map-150m.

### HASARD Specific Parameters
The following parameters are specific of the HASARD Processor.

* Blob Size / Noise Reduction ("blobsize", "NOISE_REDUCTION"): Noise Reduction is a proxy of the Blob Size Removal. Blob size removal has a range between 20 and 150. The blobs with a number of pixel < of this number will be considered noise and removed. If blobsize is set, the value is used and noise reduction is not considered. As an alternative, the user can set the Noise Reduction between 0-100%. The system will convert it in the appropriate blobsize value
* Ashman Coefficient ("ad"):  The default value of 2.6 is general, while a higher value (e.g. 2.7) can be selected to better separate the 2 distributions
* HSBA Depth ("depth"): This is the Hierarchical Split Based Approach (HSBA) Depth parameter as defined in Chini et al. (2017).

### Advanced Options
The chain works searching pre-images and post-images with the same geometry. The system divides the bounding box provided into sub-tiles of 1x1 or 2x2 or XxX degrees. 
The dimension of the grid is defined by the parameter GRIDSTEP. So, the bbox is divided in sub-squares: this is a sample of a bounding box of 3x3 degrees with a GRID step 1,1:

0,0 | 0,1 | 0,2
1,0 | 1,1 | 1,2
2,0 | 2,1 | 2,2

The final grid would be the same, for example, for a bbox of 6x6 degrees and a grid step of 2,2. 

The system generates some intermediate files:

CODE_ORBIT_YYYY-MM-DD_TILEX_TILEY.tif: mosaic of the pre-processed S1 images available of ORBIT that is located in the sub-tile X,Y
CODE_ORBIT_YYYY-MM-DD_TILEX_TILEY_[SUFFIX]: flood map obtained from S1 Images of ORBIT that is located in the sub-tile X,Y

For every day, for every orbit, the chain computes the new mosaic and searches the workspace for the last available S1 mosaic and, if present, flood in the same tile.
The final flood detection takes as input a S1 mosaic tile pre, a S1 mosaic tile post and, if available, a flood tile pre. It produces the new flood tile for the new date.

For example, assuming we have:
CODE_12_2022-01-01_0_0.tif: this is the S1 data for orbit 12 pre-processed in the sub-tile 0-0 of our bbox, for date 01/01/2022

Let's assume we have another image for orbit 12 on the 15/01/2022, for which the chain will produce
CODE_12_2022-01-15_0_0.tif

Now we have a match for a detection: CODE_12_2022-01-01_0_0.tif and CODE_12_2022-01-15_0_0.tif can be passed as input to the HASARD algorithm to produce 
CODE_12_2022-01-15_0_0_[SUFFIX]

For every day, the chain then mosaics all the flood sub-tiles to produce the daily:
CODE_2022-01-15_[SUFFIX]

But it will also keep the intermediate files of "the last valid" tiles.

Continuing with our example, we can assume there will be another S1 image for orbit 12 on 30/01/2022.
The chain will produce:
CODE_12_2022-01-30_0_0.tif
Then will search the pre-images with the same geometry and will find:
CODE_12_2022-01-15_0_0.tif
CODE_12_2022-01-15_0_0_[SUFFIX]

It will then run the HASARD algorithm with inputs: CODE_12_2022-01-15_0_0.tif, CODE_2022-01-15_[SUFFIX] and CODE_12_2022-01-30_0_0.tif to produce
CODE_12_2022-01-30_0_0_[SUFFIX]

At this poimt, the system will delete the CODE_12_2022-01-01_0_0.tif and CODE_12_2022-01-01_0_0_[SUFFIX], because they are now "old" and not anymore the last tile for a given geometry (i.e. orbit).

This is also the mechanism that can be used to re-start an archive: the system will automatically find the "last available tiles" already in the workspace.

Some advanced options are:
* MAXREFERENCETILEAGE:  max distance, in days, between a pre and post tile. Sometimes it can happen that the last valid tile become too old with respect to the new one to make a real detection. In this case, we assume that after this user-defined number of days, the reference tile is too old and should not be considered. In theory, with S1, we should have a reference tile every 15 days at maximum.
* FORCE_RE_RUN: usually the chain will re-use all the files in the workspace. However, when setting this flag to true, it forces to re-create all the existing files.
* PROVIDER: Data Provider to use in WASDI. AUTO by default.