# Plot MDT with bathymetry contours and time series of ECCO, Alt and TG
rm *.ps
ps=bathy_mdt.ps
gmt set PS_MEDIA=Custom_5.2cx4.6c
gmt set PS_PAGE_ORIENTATION=portrait

gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=thinner,20/20/20
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1p
gmt set MAP_TITLE_OFFSET=-0.18c

gmt set MAP_TICK_LENGTH_PRIMARY=0.04c
gmt set MAP_GRID_PEN_PRIMARY=thinner,lightgrey
gmt set MAP_ANNOT_OFFSET_PRIMARY=1.5p

gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,40/40/40
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,40/40/40
gmt set FONT_LABEL                     = 8p,Hind-Light,40/40/40
gmt set FONT_LOGO                      = 8p,Hind-Light,40/40/40
gmt set FONT_TITLE                     = 8p,Hind-Light,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait

gmt makecpt -CcbcOrBu.cpt -I -T-0.5/0.5/0.025 -M -D > mdt.cpt
color0=60/60/60
color1=#4e79a7
color2=#f28e2b
color3=#e15759


# mdt
gmt grdimage mdt.nc?fld -JQ0/4.5c -R260/330/-5/45 -X0.5c -Y1.0c -K -BWeSn+t'Mean dynamic topography' -By20f20 -Bx20f20 -E400 -nl -Cmdt.cpt > $ps
gmt grdcontour bathy.nc?fld -R -J -K -O -C+500 -W0.6p,50/50/50 >>$ps
gmt pscoast -R -J -Di -A5000/0/2 -G170/170/170 -K -O >> $ps
echo -e "-79.9251 32.781667" | gmt psxy -R -J -O -Sc0.07i -K -G$color3 -Wthin,50/50/50 >> $ps

gmt psscale   -Dx0.25c/-0.65c+w4.0c/0.25ch+e -R -J -Cmdt.cpt -B0.2 -O -By+l'm' -S >> $ps

gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
