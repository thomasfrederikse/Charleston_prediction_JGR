rm *.ps
ps=stats_full_alt.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14.5cx6.2c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=2p
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

# Color 0: Observations/Mean seasonal cycle
color0=60/60/60 # Mean seasonal cycle, observations (black)

# Color 1: ECCO-CCSM4
color1a=#4c78a8 # ECCO rec
color1b=#9ecae9

# Color 2: ECCO-DynP
color2a=#f58518 # Hybrid CCSM ens
color2b=#ffbf79 # Hybrid CCSM mean

# Color 3 ECCO reconstruction
color3=#54a24b # 

# Color 4: Damped persistence
color4=#e45756 # Persistence

# Color 5: CCSM ssh
color5=#72b7b2 # SLA


gmt psbasemap -K -R1/12/4.5/8.0 -JX6.25c/4c -X0.8c -Y1.9c -BWeSn+t"Root mean square error" -Bx2g2+l"Months after initialization" -By0.5g0.5 > $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.5p,$color1b -O -K rmse_ccsm_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color0  -O -K rmse_mss.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color5  -O -K rmse_ccsm_sla.txt >> $ps

gmt psxy  -R -J -N -W1.5p,$color3 -O -K rmse_ecco_obs.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color4  -O -K rmse_persistence.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a  -O -K rmse_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1a -O -K rmse_ccsm_mn.txt >> $ps
echo "1 4.5 a" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps
echo "1 6.75 RMSE (cm)" |  gmt pstext -D-0.65c/0c -R -J -N -F+f8,Hind-Light+jCM+a90 -O -K >> $ps

gmt psbasemap -O -K -R1/12/-0.2/0.7 -JX6.25c/4c -X7.25c -BWeSn+t"Anomaly correlation coefficient " -Bx2g2+l"Months after initialization" -By0.2g0.2 >> $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.5p,$color1b -O -K corr_ccsm_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color4 -O -K corr_persistence.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color5  -O -K corr_ccsm_sla.txt >> $ps

gmt psxy  -R -J -N -W1.5p,$color3 -O -K corr_ecco_obs.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a  -O -K corr_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1a -O -K corr_ccsm_mn.txt >> $ps
echo "1 -0.2 b" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps
echo "1 0.25 Correlation (-)" |  gmt pstext -D-0.65c/0c -R -J -N -F+f8,Hind-Light+jCM+a90 -O -K >> $ps

gmt psbasemap -O -K -R0/1/-1.5/1.5 -JX9.5c/1.1c -X-5.25c -Y-1.9c -T100/100/1 >> $ps

echo -e "0 1 \n 0.03 1" | gmt psxy -R -J -O -K -W1.5p,$color1a >> $ps
echo -e "0.03 1 \n 0.06 1" | gmt psxy -R -J -O -K -W0.75p,$color1b >> $ps
echo "0.065 1 Hybrid prediction (ECCO-CCSM4)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0.0 0 \n 0.06 0" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.065 0 Mean seasonal cycle" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0.50 1 \n 0.56 1" | gmt psxy -R -J -O -K -W1.5p,$color2a >> $ps
echo "0.565 1 Hybrid prediction (ECCO-DynP)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0.50 0 \n 0.56 0" | gmt psxy -R -J -O -K -W1.5p,$color4 >> $ps
echo "0.565 0 Damped persistence" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0 -1 \n 0.06 -1" | gmt psxy -R -J -O -K -W1.5p,$color3 >> $ps
echo "0.065 -1 ECCO reconstruction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "0.50 -1 \n 0.56 -1" | gmt psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "0.565 -1 CCSM4 SLA prediction" | gmt pstext -R -J -F+f8+jLM -O >> $ps

gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
