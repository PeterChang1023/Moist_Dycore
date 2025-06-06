#!/bin/sh

# Define pressure levels for vertical interpolation (in hPa)
# levels="10,20,30,50,70,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,600,650,700,750,775,800,850,900,950,1000"
levels="18.5,44.64,70.85,97.06,123.28,149.49,175.70,201.91,228.13,254.34,280.55,306.76,332.97,359.19,385.40,411.61,437.82,464.04,490.25,516.46,542.67,568.89,595.10,621.31,647.52,673.74,699.95,726.16,752.37,778.59,804.80,831.01,857.22,883.44,909.65,935.86,962.08"



# Define horizontal target grid (2.5°x2.5°, from lon=-180 to 180 and lat=90 to -90)
cat > mygrid << EOF
gridtype = lonlat
xsize    = 144
ysize    = 73
xfirst   = -180
xinc     = 2.5
yfirst   = 90
yinc     = -2.5
EOF

# Loop over variables
for var in u v t
do
  echo "Processing variable: $var"
 
  infile="Dycore_500_20000day_6hourly_${var}.nc"
  hfile="Dycore_500_20000day_6hourly_${var}_hinterp.nc"
  outfile="Dycore_500_20000day_6hourly_${var}_regrid.nc"
  outfile_sorted="Dycore_500_20000day_6hourly_${var}_regrid_sorted.nc"
  
  # rm $outfile 
  # 1. Horizontal interpolation
  cdo remapbil,mygrid $infile $hfile

  # 2. Vertical interpolation
  cdo intlevel,$levels $hfile $outfile

  # 3. Reset time axis (6-hourly from 1979-01-01)
  cdo settaxis,1979-01-01,00:00:00,6hour $outfile ${outfile%.nc}_taxis.nc
  mv ${outfile%.nc}_taxis.nc $outfile

  # 4. Force longitude range [-180, 180) and ensure lat is descending
  cdo invertlat $outfile tmp1.nc
  cdo sellonlatbox,-180,180,-90,90 tmp1.nc $outfile_sorted
  mv $outfile_sorted $outfile
  rm -f tmp1.nc $hfile
done

# Clean up
rm -f mygrid

echo "✅  All variables processed, regridded, and ordered lon[-180~180], lat[90~-90]."

