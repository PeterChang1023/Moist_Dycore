#!/bin/sh

# Create a target grid file with correct names and attributes
cat > mygrid << EOF
gridtype = lonlat
gridsize = 10512
xsize = 144
ysize = 73
xname = longitude
xlongname = "longitude"
xunits = "degrees_east"
yname = latitude
ylongname = "latitude"
yunits = "degrees_north"
xfirst = -180
xinc = 2.5
yfirst = 90
yinc = -2.5
EOF

# Loop over variables
for var in u v t
do
  infile="Dycore_500_20000day_6hourly_${var}_regrid.nc"
  outfile="Dycore_500_20000day_6hourly_${var}_final.nc"

  echo "✅ Regridding $var to target grid with correct dimension/attribute names..."
  cdo -f nc4 remapcon,mygrid $infile $outfile
done

# Clean up
rm -f mygrid

echo "✅ Done: latitude/longitude renamed and standardized using CDO."
