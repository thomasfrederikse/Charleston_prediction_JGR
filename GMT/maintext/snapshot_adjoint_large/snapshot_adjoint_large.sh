rm *.ps
ps=snapshot_adjoint_large.ps
gmt set PS_MEDIA=Custom_16.6cx11.25c

gmt set MAP_LABEL_OFFSET=1p
gmt set MAP_ANNOT_OFFSET=1p

gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_TICK_LENGTH_PRIMARY=0.06c
gmt set MAP_GRID_PEN_PRIMARY=default,lightgrey
gmt set MAP_TITLE_OFFSET=-0.2c
gmt set FORMAT_GEO_MAP=D
gmt set MAP_FRAME_PEN=0.5p,40/40/40
gmt set FONT_ANNOT_PRIMARY             = 7p,Hind-Light,black
gmt set FONT_ANNOT_SECONDARY           = 7p,Hind-Light,black
gmt set FONT_LABEL                     = 7p,Hind-Light,black
gmt set FONT_LOGO                      = 7p,Hind-Light,black
gmt set FONT_TITLE                     = 7p,Hind-Light,black
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_COLOR_MODEL=RGB
gmt set MAP_FRAME_TYPE=plain

J1=Q3.5c
R1=260/360/0/60
color1=100/100/100
xjump=3.6c
yjump=-2.8c
xback=-10.8c

fn_sens=/Users/tfrederi/Data/ECCO/RISE/results/Charleston/sens_jan_interp.nc
gmt makecpt -CcbcOrBu.cpt -I -T-0.1/0.1/0.005 -M -D > empmr.cpt
gmt makecpt -CcbcPRGn.cpt -I -T-4e-7/4e-7/2e-8 -M -D > qnet.cpt
#gmt makecpt -CcbcPuOr.cpt -I -T-1e-3/1e-3/5e-5 -M -D > tauu.cpt
gmt makecpt -CcbcRdBu.cpt -I -T-4e-4/4e-4/2e-5 -M -D > tauu.cpt
gmt makecpt -CcbcRdBu.cpt -I -T-4e-4/4e-4/2e-5 -M -D > tauv.cpt

###
# Zonal wind
gmt grdimage $fn_sens?tauu[12] -J$J1 -R$R1 -X0.8c -Y8.75c -K -BWrst+t"0 months" -By20f20 -Bx30f30 -Ctauu.cpt > $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
echo "260 30 Zonal wind stress" |  gmt pstext -D-0.7c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps
gmt grdimage $fn_sens?tauu[8] -J$J1 -R$R1 -X$xjump -O -K -Bwrst+t"4 months" -By20f20 -Bx30f30 -Ctauu.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?tauu[4] -J$J1 -R$R1 -X$xjump -O -K -Bwrst+t"8 months" -By20f20 -Bx30f30 -Ctauu.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?tauu[0] -J$J1 -R$R1 -X$xjump -O -K -Bwrst+t"12 months" -By20f20 -Bx30f30 -Ctauu.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt psscale -Dx3.85c/0.0c+w2.1c/0.25+e -R -J -Np -Ctauu.cpt -B2e-4 -By+l'm N@+-1@+ m@+2@+' -O -K -S --FORMAT_FLOAT_MAP=%.0e >> $ps

# Meridional wind
gmt grdimage $fn_sens?tauv[12] -J$J1 -R$R1 -X$xback -Y$yjump -O -K -BWrst -By20f20 -Bx30f30 -Ctauv.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
echo "260 30 Meridional wind stress" |  gmt pstext -D-0.7c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps
gmt grdimage $fn_sens?tauv[8] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Ctauv.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?tauv[4] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Ctauv.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?tauv[0] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Ctauv.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt psscale -Dx3.85c/0.0c+w2.1c/0.25+e -R -J -Np -Ctauv.cpt -B2e-4 -By+l'm N@+-1@+ m@+2@+' -O -K -S --FORMAT_FLOAT_MAP=%.0e >> $ps

# Heat
gmt grdimage $fn_sens?qnet[12] -J$J1 -R$R1 -X$xback -Y$yjump -O -K -BWrst -By20f20 -Bx30f30 -Cqnet.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
echo "260 30 Heat flux" |  gmt pstext -D-0.7c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps
gmt grdimage $fn_sens?qnet[8] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Cqnet.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?qnet[4] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Cqnet.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?qnet[0] -J$J1 -R$R1 -X$xjump -O -K -Bwrst -By20f20 -Bx30f30 -Cqnet.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt psscale -Dx3.85c/0.0c+w2.1c/0.25+e -R -J -Np -Cqnet.cpt -B2e-7 -By+l'm W@+-1@+ m@+2@+' -O -K -S --FORMAT_FLOAT_MAP=%.0e >> $ps

# Freshwater
gmt grdimage $fn_sens?empmr[12] -J$J1 -R$R1 -X$xback -Y$yjump -O -K -BWrSt -By20f20 -Bx30f30 -Cempmr.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
echo "260 30 Freshwater flux" |  gmt pstext -D-0.7c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps
gmt grdimage $fn_sens?empmr[8] -J$J1 -R$R1 -X$xjump -O -K -BwrSt -By20f20 -Bx30f30 -Cempmr.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?empmr[4] -J$J1 -R$R1 -X$xjump -O -K -BwrSt -By20f20 -Bx30f30 -Cempmr.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt grdimage $fn_sens?empmr[0] -J$J1 -R$R1 -X$xjump -O -K -BwrSt -By20f20 -Bx30f30 -Cempmr.cpt >> $ps
gmt pscoast -R -J -Dl -A1000/0/2 -G170/170/170  -O -K >> $ps
echo "-79.925 32.78" | gmt psxy -R -J -O -Sc0.06i -K -G$color1 -Wthin >> $ps
gmt psscale -Dx3.85c/0.0c+w2.1c/0.25+e -R -J -Np -Cempmr.cpt -B0.05 -By+l'm kg@+-1@+ m@+2@+ s' -O -S >> $ps

gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps


