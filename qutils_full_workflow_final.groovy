import qupath.lib.gui.measure.ObservableMeasurementTableData

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
double pixelWidth = cal.getPixelWidthMicrons()
double pixelHeight = cal.getPixelHeightMicrons()

// !!Normalize to project-specific vectors!!

def imgName = getCurrentServer().getMetadata().getName()

double amtTotal = 0 //Total number of detections
double amtImmune = 0 //Total number of lymphocyte-classified detections
double allCellsInSpotTotal = 0 //Total number of detections in identified TIL hotspot
double allCellsInSpotMarg = 0 //Total number of detections in identified TIL hotspot that falls within calculated invasive margin area
double immuneInSpotTotal = 0 //Total number of detections classed as lymphocyte in TIL hotspot (subset of allCellsInSpotTotal)
double immuneInSpotMarg = 0 //Total number of detections classed as lymphocyte in TIL hotspot that falls within calculated invasive margin area (subset of allCellsInSpotMarg)


def avgImmuneTotal = [] //List for calculating average number of detected lymphocytes across all TIL hotspots
def avgImmuneMarg = [] //List for calculating average number of detected lymphocytes across all TIL hotspots overlapping estimated invasive margin area

double marg_tils = 0 //Total number of detections classed as lymphocytes in the estimated invasive margin area
double marg_tum = 0 //Total number of detections classed as tumor in the estimated invasive margin area
double marg_stroma = 0 //Total number of detections classed as stroma in the estimated invasive margin area
//Downstream calculation for marginal TILs: marg_tils/(marg_tils + marg_stroma); estimated area often results in poor segmentation of the exact boundary lines for invasive margin compartments, most consistent calculation approaching manual estimations of stromal TILs


// Don't use this if image annotations are (ideally) pre-existing
//createFullImageAnnotation(true)
//getSelectedObjects().each {
//    it.setPathClass(getPathClass("Positive"))
//}


   
   
//Breaks tissue annotation into tiles, puts them in callable list:

selectAnnotations()
getSelectedObjects().each {
    it.setPathClass(getPathClass("Positive")) //Can change to whatever class desired, just needs to be consistent so the script can refer to it by classification
    }
mergeSelectedAnnotations()
runPlugin('qupath.lib.algorithms.TilerPlugin', '{"tileSizeMicrons":1000.0,"trimToROI":true,"makeAnnotations":true,"removeParentAnnotation":false}')

def allAnnotations = getAnnotationObjects()
def annotations = []
allAnnotations.each { annotation ->
    if (annotation.getName()?.startsWith('Tile ')) {
        annotations.add(annotation)}
    }

//Detect cells in each tile - useful for avoiding possible memory errors by detecting across entire tissue ROI

for(anno in annotations) {
    selectObjects { p -> p == anno }
        runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield":"Hematoxylin OD","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.5,"minAreaMicrons":10.0,"maxAreaMicrons":400.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"cellExpansionMicrons":5.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')
    //Dump cache memory, may help avoid crashing in particularly large images
    javafx.application.Platform.runLater {
    getCurrentViewer().getImageRegionStore().cache.clear()
    System.gc()
    }
}

//Add QuTILs required classification features

selectCells();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":1.0,"region":"NUCLEUS","tileSizeMicrons":10.0,"colorOD":false,"colorStain1":true,"colorStain2":true,"colorStain3":false,"colorRed":false,"colorGreen":fase,"colorBlue":false,"colorHue":false,"colorSaturation":false,"colorBrightness":false,"doMean":true,"doStdDev":true,"doMinMax":false,"doMedian":false,"doHaralick":true,"haralickDistance":1,"haralickBins":32}')
selectObjectsByClassification("Positive")
runPlugin('qupath.lib.plugins.objects.SmoothFeaturesPlugin', '{"fwhmMicrons":20.0,"smoothWithinClasses":false}')
runPlugin('qupath.lib.plugins.objects.SmoothFeaturesPlugin', '{"fwhmMicrons":40.0,"smoothWithinClasses":false}')
      
runObjectClassifier("tils-final"); //Classifies cells by QuTILs categories: Tumor, Stroma, Immune, or 'Other'
    
javafx.application.Platform.runLater {
getCurrentViewer().getImageRegionStore().cache.clear()
System.gc()
}

//First find immune-related hotspots

selectObjectsByClassification("Positive")
findDensityMapHotspots("tils-map", 0, 50000, 10.000000, false, false) //Generates estimated TILs hotspots - very rough, may not be useful

def hotspotData_current = getAnnotationObjects().findAll{it.getName()?.startsWith('Hotspot')}
def removeBad = [] //List really poor identified hotspots for deletion

def imageData_tils = getCurrentImageData()
def ob1 = new ObservableMeasurementTableData()
ob1.setImageData(imageData_tils, hotspotData_current)

//Ensure TILs hotspots actually have required number of TILs

hotspotData_current.each {immuneSpots ->
    if (ob1.getNumericValue(immuneSpots, 'Num Immune cells') < 10) {
       removeBad << immuneSpots 
    }
}

removeObjects(removeBad, false)


def hotspotData_new = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Immune cells: Immune cells")}

def imageData_new = getCurrentImageData()
def ob2 = new ObservableMeasurementTableData()
ob2.setImageData(imageData_new, hotspotData_new)

hotspotData_new.each {immuneSpotsFinal ->
    allCellsInSpotTotal += ob2.getNumericValue(immuneSpotsFinal, 'Num Detections')
    immuneInSpotTotal += ob2.getNumericValue(immuneSpotsFinal, 'Num Immune cells')
    avgImmuneTotal << ob2.getNumericValue(immuneSpotsFinal, 'Num Immune cells')
    }


//Get total whole slide counts:

def totalCells = getDetectionObjects().findAll()
amtTotal = totalCells.size()

def immuneCount = getDetectionObjects().findAll {it.getPathClass() == getPathClass("Immune cells")}
amtImmune = immuneCount.size()

//Total hotspots and average TILs count across total detected hotspots

double amtHotSpotsTotal = hotspotData_new.size()

//get average TILs in hotspots else 0 if no hotspots

def getSpotsAveTotal = avgImmuneTotal ? (avgImmuneTotal.sum() / avgImmuneTotal.size()) : 0

//Coarsely estimate invasive margin by finding dense tumor pockets
//Join them, expand outward by 500µm then dilate inward by 1000 µm

createAnnotationsFromDensityMap("tumor-map", [0: 1000.0], "Tumor cells: Tumor cells", "SPLIT")

tumorAnno = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells")}

if (tumorAnno.isEmpty()) {
    createAnnotationsFromDensityMap("tumor-map", [0: 500.0], "Tumor cells: Tumor cells", "SPLIT")
}

tumorAnno = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells")}

if (tumorAnno.isEmpty()) {
        createAnnotationsFromDensityMap("tumor-map", [0: 100.0], "Tumor cells: Tumor cells", "SPLIT")
}

tumorAnno = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells")}

def ob3 = new ObservableMeasurementTableData();
    ob3.setImageData(getCurrentImageData(), tumorAnno);
    
//Get area to sort by size

def measureList = []

tumorAnno.each{tA -> 
    annoArea = ob3.getNumericValue(tA, "Area µm^2")
    measureList << annoArea
    }

measureList.sort {b, a -> a <=> b}
def tumorThreshold = 0

//If biggest tumor section is large (300k), keep, otherwise keep top 5 tumor:

if (measureList[0] >= 300000) {
    tumorThreshold = measureList[1]
} else if (measureList.size() >= 5){
    tumorThreshold = measureList[4]
} else {
    tumorThreshold = measureList[-1]
}

def removeSmall = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells") && it.getROI().getScaledArea(pixelWidth, pixelHeight) < tumorThreshold}
removeObjects(removeSmall, true)

//Expand outward and dilate inward leaving simulated invasive margin

selectObjectsByClassification("Tumor cells: Tumor cells")
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons":500.0,"lineCap":"ROUND","removeInterior":false,"constrainToParent":true}')
selectObjectsByClassification("Tumor cells: Tumor cells")
mergeSelectedAnnotations()
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons":-1000.0,"lineCap":"ROUND","removeInterior":true,"constrainToParent":true}')
getSelectedObjects().each {
    it.setPathClass(getPathClass("AllTumor"))
}
resetSelection()


def margArea = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells")}
    def ob4 = new ObservableMeasurementTableData();
    ob4.setImageData(getCurrentImageData(), margArea);


margArea.each { cells -> 
    marg_tils += ob4.getNumericValue(cells, 'Num Immune cells')
    marg_tum += ob4.getNumericValue(cells, 'Num Tumor cells')
    marg_stroma += ob4.getNumericValue(cells, 'Num Stroma cells')
    }
    


invMargin = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Tumor cells: Tumor cells")}
annoSpots = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Immune cells: Immune cells")}

double margSpots = 0
def nonInvSpots = []

// !!Should be only one merged tumor anno!! -- else will remove all hotspots

invMargin.each {invMarg -> 
    def marginalROI = invMarg.getROI().getGeometry()
    
    annoSpots.each{immSpot -> 
        def spotROI = immSpot.getROI().getGeometry()
        def intersectTest = spotROI.intersection(marginalROI)
        //if intersect w/ invasive margin, keep; if not, remove
        if (!intersectTest.isEmpty()) {
            margSpots++
        }
        if (intersectTest.isEmpty()) {
            nonInvSpots << immSpot
        }

    }
}

removeObjects(nonInvSpots)

//Repeat TILs hotspots feature calculations but only for those in range of simulated invasive margin

def hotspotData_margspot = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Immune cells: Immune cells")}

def imageData_margspot = getCurrentImageData()
def ob5 = new ObservableMeasurementTableData()
ob5.setImageData(imageData_margspot, hotspotData_margspot)

if (hotspotData_margspot.size() > 0) {
    hotspotData_margspot.each {immuneSpotsMarginal ->
        allCellsInSpotMarg += ob5.getNumericValue(immuneSpotsMarginal, 'Num Detections')
        immuneInSpotMarg += ob5.getNumericValue(immuneSpotsMarginal, 'Num Immune cells')
        avgImmuneMarg << ob5.getNumericValue(immuneSpotsMarginal, 'Num Immune cells')
    }
}

double amtHotSpotsMarg = hotspotData_margspot.size()

def getSpotsAveMarg = avgImmuneMarg ? (avgImmuneMarg.sum() / avgImmuneMarg.size()) : 0

double est_stromal_invasive_marginal_tils = marg_tils / (marg_tils + marg_stroma) * 100

def resultsPath = buildFilePath(PROJECT_BASE_DIR, "tils_results.txt")
new File(resultsPath).withWriterAppend {writer ->
    writer.writeLine("$imgName    $amtTotal    $amtImmune    $marg_tils    $marg_tum    $marg_stroma    $amtHotSpotsTotal    $allCellsInSpotTotal    $immuneInSpotTotal    $getSpotsAveTotal    $amtHotSpotsMarg    $allCellsInSpotMarg    $immuneInSpotMarg    $getSpotsAveMarg    $est_stromal_invasive_marginal_tils")
}

print(imgName + " - Job's done.")