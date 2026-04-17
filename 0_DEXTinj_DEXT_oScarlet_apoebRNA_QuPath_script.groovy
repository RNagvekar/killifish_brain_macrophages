//This script will detect all DAPI+ cells.

////////////////////////////////////////////////////////////////////
//Parameters:
////////////////////////////////////////////////////////////////////


//Rename channels to be informative:
setChannelNames(
     'DAPI',
     'DEXT',
     'oScarlet',
     'apoebRNA'
)

/*Set channel colors
def colors = [
    getColorRGB(255,0,0),
    getColorRGB(0,255,0),
    getColorRGB(0,0,255),
    getColorRGB(255,255,255)
    ]
setChannelColors( *colors )
*/

//Set brightness:
setChannelDisplayRange(0,0,20000)
setChannelDisplayRange(1,0,20000)
setChannelDisplayRange(2,0,20000)
setChannelDisplayRange(3,0,20000)

//Script to detect cells by DAPI

resetSelection();
selectDetections();
clearSelectedObjects();
clearAnnotations();
  
//Select all annotations
//selectAnnotations();
createFullImageAnnotation(true)

//Nuclei detection
    //sigma: 1.5 um
    //expand: 2.0 um
    //min area: 10 um 
    //max area: 400 um
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 2,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 500,  "watershedPostProcess": true,  "cellExpansionMicrons": 2.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

fireHierarchyUpdate()

//Access measurements for entire project:
//Measure > Export measurements
// select "Annotations" for counts
// select "Detections" for cell level
// select "Cells" for signal and cell area
//Save.