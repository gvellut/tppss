function DemToTIFF([Parameter(Mandatory=$true)][string]$dataDirPath, [Parameter(Mandatory=$true)][string]$name, [Parameter(Mandatory=$true)][string]$outputDirPath) {
    $tiffListFileName = "tiffs.txt"
    New-Item -Path "$dataDirPath" -Name "$tiffListFileName" -ItemType "file" -Force

    Get-ChildItem $dataDirPath -Filter *.asc | 
    Foreach-Object {
        $outputPath = $_.DirectoryName + "\" + $_.BaseName + ".tif"
        gdal_translate -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" "$($_.FullName)" "$outputPath"
        Add-Content "$dataDirPath\$tiffListFileName" "$outputPath"
    }

    gdalbuildvrt "$dataDirPath\mosaic.vrt" -input_file_list "$dataDirPath\$tiffListFileName"
    gdalwarp -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" -overwrite "$dataDirPath\mosaic.vrt" "$outputDirPath\$name.tif"  -t_srs EPSG:5698
}

function MergeTIFF([Parameter(Mandatory=$true)][string]$dataDirPath) {
    $tiff_paths = (ls $dataDirPath *.tif).FullName

    $rgf93_tiff_path = "C:\Users\gvellut\Documents\projects\sunlight\dem_rgf93.tif"
    $wgs84_tiff_path = "C:\Users\gvellut\Documents\projects\sunlight\dem_wgs84_b.tif"

    python C:\Users\gvellut\anaconda3\envs\calcsoleil\Scripts\gdal_merge.py -o  "$rgf93_tiff_path" -co "COMPRESS=LZW" -co "TILED=YES" $tiff_paths

    # To transform a DEM from geoid elevations (using RGF93 / Lambert-93 + NGF-IGN69) to WGS84 ellipsoidal heights
    gdalwarp -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" -overwrite "$rgf93_tiff_path" "$wgs84_tiff_path"  -t_srs EPSG:4979
}