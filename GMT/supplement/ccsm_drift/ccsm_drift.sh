rm *.ps
ps=ccsm_drift.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14.5cx10.0c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1p
gmt set MAP_TITLE_OFFSET=0.05c
gmt set MAP_TICK_LENGTH_PRIMARY=0.03c
gmt set MAP_TICK_PEN=thinnest,40,40,40
gmt set MAP_GRID_PEN_PRIMARY=thinnest,lightgrey
gmt set MAP_ANNOT_OFFSET_PRIMARY=1.0p
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,black
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,black
gmt set FONT_LABEL                     = 8p,Hind-Light,black
gmt set FONT_LOGO                      = 8p,Hind-Light,black
gmt set FONT_TITLE                     = 8p,Hind-Light,black

# COLORS
color0=80/80/80
color1=#4c78a8
color2=#f58518
gmt makecpt -CcbcOrBu.cpt -I -T-100/100/5 -M -D > drift.cpt

# Panel a: drift in Charleston
gmt psbasemap -K -R1995/2017/-40/60 -JX7.75c/3.5c -X0.8c -Y6.2c -BWeSn -Bx4g4 -By20g20+l"Charleston sea level (cm)" > $ps
gmt psxy  -R -J -N -W1.0p,$color0 -O -K rsl_obs.txt >> $ps
for i in {1995..2016}
do
   gmt psxy  -R -J -N -W1.25p,$color1 -O -K ${i}.txt >> $ps
   gmt psxy  -R -J -N -W1.25p,$color2 -O -K ${i}_noseas.txt >> $ps
done
gmt psxy  -R -J -N -Sc0.15c -G$color1 -O -K startpoint.txt >> $ps
gmt psxy  -R -J -N -Sh0.15c -G$color2 -O -K noseas_startpoint.txt >> $ps
echo "1995 -40 a" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

# Panel b: drift map
gmt psbasemap -O -K -R0/0.9/-1.5/1.5 -JX6.75c/1.0c -X0.5c -Y-1.3c -T100/100/1  >> $ps
echo -e "0.03 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color1 >> $ps
echo -e "0.03 1 \n 0.1 1" | gmt psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.11 1 CCSM4 SSH prediction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "0.03 0" | gmt psxy -R -J -O -K -Sh0.15c -G$color2 >> $ps
echo -e "0.03 0 \n 0.1 0" | gmt psxy -R -J -O -K -W1.5p,$color2 >> $ps
echo "0.11 0 CCSM4 SSH prediction (Seasonal cycle removed)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
 echo -e "0.03 -1 \n 0.1 -1" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.11 -1 Observed sea level (Seasonal cycle removed)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

gmt grdimage ccsm_drift_map.nc?ydiff -Rd -JX5.0c/3.5c -X7.95c -Y1.3c -BWeSn -By40f40 -Bx60f60 -E400 -nl -Cdrift.cpt -O -K >> $ps
gmt pscoast -R -J -Dc -A5000/0/1 -G170/170/170 -K -O >> $ps
gmt psscale -Dx0.5c/-0.8c+w4.0c/0.25ch+e -R -J -Cdrift.cpt -B30 -O -K -By+l'cm' -S >> $ps
echo "-180 -90 b" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -t50 -O -K >> $ps
echo "-180 -90 b" | gmt pstext -D0.125c/0.125c -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

# Panel C: RMSE
gmt set MAP_LABEL_OFFSET=3p
gmt psbasemap -O -K -R1/12/6.4/8.6 -JX6c/3.5c -X-8.45c -Y-5.5c -BWeSn -Bx2g2+l"Prediction lead time (months)" -By0.4g0.4+l"Root Mean Square Error (cm)" >> $ps
gmt psxy  -R -J -N -W1.5p,$color0  -O -K rmse_mss.txt >> $ps
gmt psxy  -R -J -N -Sc0.2c -G$color0  -O -K rmse_mss.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1  -O -K rmse_sla.txt >> $ps
gmt psxy  -R -J -N -Sh0.2c -G$color1  -O -K rmse_sla.txt >> $ps
echo -e "6.0 8.45" | gmt psxy -R -J -O -K -Sh0.2c -G$color1 >> $ps
echo -e "5.5 8.45 \n 6.5 8.45" | gmt psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "6.55 8.45 CCSM4 SSH prediction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "6.0 8.25" | gmt psxy -R -J -O -K -Sc0.2c -G$color0 >> $ps
echo -e "5.5 8.25 \n 6.5 8.25" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "6.55 8.25 Mean seasonal cycle" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps


echo "1 6.4 c" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1/12/-0.1/0.3 -JX6c/3.5c -X7.5c -BWeSn -Bx2g2+l"Prediction lead time (months)" -By0.1g0.1+l"Anomaly correlation coefficient (-)" >> $ps
gmt psxy  -R -J -N -W1.5p,$color1  -O -K corr_sla.txt >> $ps
gmt psxy  -R -J -N -Sh0.2c -G$color1  -O -K corr_sla.txt >> $ps
echo "1 -0.1 d" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O >> $ps


gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

