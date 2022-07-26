rm *.ps
ps=example_rec_tseries_tg.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14.5cx10.9c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1.5p
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

color0=80/80/80 # Observed 
color1=#4c78a8 # ECCO-CCSM4
color2=#f58518 # ECCO-DynP
# color3=#54a24b # 
# color4=#e45756 # 
color5=#72b7b2 # CCSM SLA

gmt psbasemap -K -R1995/2017/-17.5/18.5 -JX13.2c/4.0c -X0.8c -Y6.4c -BWeSn+t"Prediction one month ahead" -Bx3g3 -By5g5+l"Sea level (cm)" > $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.5p,$color5 -O -K sla_pred_ens_${i}_1m.txt -t60 >> $ps
   gmt psxy  -R -J -N -W0.5p,$color1 -O -K obs_pred_ens_${i}_1m.txt -t60 >> $ps
done
gmt psxy  -R -J -N -W1.25p,$color5 -O -K sla_pred_1m.txt >> $ps
gmt psxy  -R -J -N -W1.25p,$color0 -O -K obs_1m.txt >> $ps
gmt psxy  -R -J -N -W1.25p,$color2 -O -K obs_seas_1m.txt >> $ps
gmt psxy  -R -J -N -W1.25p,$color1 -O -K obs_pred_1m.txt >> $ps
echo "1995 -17.5 a" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1995/2017/-17.5/18.5 -JX13.2c/4.0c  -Y-5.0c -BWeSn+t"Prediction four months ahead" -Bx3g3+l"Time (years)" -By5g5+l"Sea level (cm)" >> $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.5p,$color5 -O -K sla_pred_ens_${i}_4m.txt -t60 >> $ps
   gmt psxy  -R -J -N -W0.5p,$color1 -O -K obs_pred_ens_${i}_4m.txt -t60 >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color5 -O -K sla_pred_4m.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K obs_4m.txt >> $ps
gmt psxy  -R -J -N -W1.25p,$color2 -O -K obs_seas_4m.txt >> $ps
gmt psxy  -R -J -N -W1.25p,$color1 -O -K obs_pred_4m.txt >> $ps
echo "1995 -17.5 b" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps


gmt psbasemap -O -K -R0/2/0.5/2.5 -JX9.2c/0.75c -X2c -Y-1.4c -T100/100/1  >> $ps
echo -e "0 2 \n 0.1 2" | gmt psxy -R -J -O -K -W1.25p,$color0 >> $ps
echo "0.105 2 Observed SLA (Tide gauge)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0.0 1 \n 0.05 1" | gmt psxy -R -J -O -K -W1.25p,$color5 >> $ps
echo -e "0.05 1 \n 0.1 1" | gmt psxy -R -J -O -K -W0.5p,$color5 >> $ps
echo "0.105 1 CCSM4 SLA prediction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "1 2 \n 1.05 2" | gmt psxy -R -J -O -K -W1.25p,$color1 >> $ps
echo -e "1.05 2 \n 1.1 2" | gmt psxy -R -J -O -K -W0.5p,$color1 >> $ps
echo "1.105 2 Hybrid prediction (ECCO-CCSM4)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "1 1 \n 1.10 1" | gmt psxy -R -J -O -K -W1.25p,$color2 >> $ps
echo "1.105 1 Hybrid prediction (ECCO-DynP)" | gmt pstext -R -J -F+f8+jLM -O >> $ps


gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

