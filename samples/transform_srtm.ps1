function DemToTIFF([Parameter(Mandatory=$true)][string]$dataPath, [Parameter(Mandatory=$true)][string]$outputDirPath) {
    $dataPathObj = Get-Item $dataPath
    $outputPath = $outputDirPath + "\" + $dataPathObj.BaseName + ".tif"
    gdal_translate -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" "$dataPath" "$outputPath"
    
}

function MergeTIFF([Parameter(Mandatory=$true)][string]$dataDirPath) {
    $tiff_paths = (ls $dataDirPath *.tif).FullName
    #$tiff_paths = '"{0}"' -f ($tiff_paths -Join '" "')

    $wgs84_tiff_path = "C:\Users\gvellut\Documents\projects\sunlight\dem_wgs84_srtm.tif"

    python C:\Users\gvellut\anaconda3\envs\calcsoleil\Scripts\gdal_merge.py -o  "$wgs84_tiff_path" -co "COMPRESS=LZW" -co "TILED=YES" $tiff_paths

}