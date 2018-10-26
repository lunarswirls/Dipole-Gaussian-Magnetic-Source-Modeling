#!/bin/bash
in2="ReinerGamma_z_i0d-20.in"
out2="ReinerGamma_z_i0d-20_e30.xyz"
gfortran -O -o moonsphere5_gauss2_32_brtp moonsphere5_gauss2_32_brtp.f
cp gauss2_32.data fort.1
#
cp $in2 fort.4
./moonsphere5_gauss2_32_brtp
mv fort.2 $out2
#
rm fort.*
rm sphereiiscratch
#
in="ReinerGamma_z_i0d-20_e30"
nx=151
ny=151
xmin=295.0
xmax=310.0
xinc=0.1
ymin=0.0
ymax=15.0
yinc=0.1
#
#		@(#)map.bash	1.1  09/21/98
#
# Purpose:	Plot gridded data
# GMT progs:	gmtset grd2cpt grdgradient grdimage makecpt psscale pstext
# Unix progs:	cat rm
#
#gmtset HEADER_FONT_SIXE 30 OBLIQUE_ANOTATION 0 DEGREE_FORMAT 0
gmtset PROJ_LENGTH_UNIT inch PS_MEDIA letter
xyz2grd $in.xyz -G$in.grd -R$xmin/$xmax/$ymin/$ymax -I$xinc
#grd2cpt $in.grd > $in.cpt 
grdimage $in.grd -K -R$xmin/$xmax/$ymin/$ymax -JL307.5/-89/0/15/4i -Ba5g1f1NEsw -V -C$in.cpt -X2 -Y5 -P > $in.eps 
psscale -C$in.cpt -D+5.0/-1.5/10.00/0.5h -L -O -V  -B:nT: >> $in.eps 
\rm -f .gmtcommands
