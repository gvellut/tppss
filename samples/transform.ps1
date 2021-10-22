param([string]$dataDirPath, [string]$name, [string]$outputDirPath)

$tiffListFileName = "tiffs.txt"
New-Item -Path "$dataDirPath" -Name "$tiffListFileName" -ItemType "file" -Force

Get-ChildItem $dataDirPath -Filter *.asc | 
Foreach-Object {
    $outputPath = $_.DirectoryName + "\" + $_.BaseName + ".tif"
    gdal_translate -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" "$($_.FullName)" "$outputPath"
    Add-Content "$dataDirPath\$tiffListFileName" "$outputPath"
}

gdalbuildvrt "$dataDirPath\mosaic.vrt" -input_file_list "$dataDirPath\$tiffListFileName"

# To transform a DEM from geoid elevations (using RGF93 / Lambert-93 + NGF-IGN69) to WGS84 ellipsoidal heights
gdalwarp -of GTiff -co "COMPRESS=LZW" -co "TILED=YES" -overwrite "$dataDirPath\mosaic.vrt" "$outputDirPath\$name.tif"  -s_srs EPSG:5698 -t_srs EPSG:4979

# TODO mosaic 73 + 74 + 38 + 01 (+ geneve)
