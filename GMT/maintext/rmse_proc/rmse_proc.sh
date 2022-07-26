rm *.ps
ps=rmse_proc.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14cx9.5c
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


gmt psbasemap -K -R1/12/0/0.5 -JX6.2c/4.0c -X0.7c -Y5.1c -BWesn+t"Freshwater flux" -Bx2g2 -By0.1g0.1+l"Root mean square error (cm)" > $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.75p,$color1b -O -K rmse_empmr_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color1a -O -K rmse_empmr_mean.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a -O -K rmse_empmr_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K rmse_empmr_seas.txt >> $ps
echo -e "3.2 0.11 \n 3.7 0.11" | gmt psxy -R -J -O -K -W1.5p,$color1a >> $ps
echo -e "3.7 0.11 \n 4.2 0.11" | gmt psxy -R -J -O -K -W0.75p,$color1b >> $ps
echo "4.3 0.11 Hybrid prediction (ECCO-CCSM4)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "3.2 0.07 \n 4.2 0.07" | gmt psxy -R -J -O -K -W1.5p,$color2a >> $ps
echo "4.3 0.07 Hybrid prediction (ECCO-DynP)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "3.2 0.03 \n 4.2 0.03" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "4.3 0.03 Mean seasonal cycle" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo "1 0 a" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1/12/0.4/2.4 -JX6.2c/4.0c -X6.9c -BWesn+t"Heat flux" -Bx2g2 -By0.4g0.4 >> $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.75p,$color1b -O -K rmse_qnet_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color1a -O -K rmse_qnet_mean.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a -O -K rmse_qnet_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K rmse_qnet_seas.txt >> $ps
echo "1 0.4 b" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1/12/2/4 -JX6.2c/4.0c -X-6.9c -Y-4.5c -BWeSn+t"Zonal wind" -Bx2g2+l"Months after initialization" -By0.4g0.4+l"Root mean square error (cm)" >> $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.75p,$color1b -O -K rmse_tauu_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color1a -O -K rmse_tauu_mean.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a -O -K rmse_tauu_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K rmse_tauu_seas.txt >> $ps
echo "1 2 c" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1/12/1.4/2.4 -JX6.2c/4.0c -X6.9c -BWeSn+t"Meridional wind" -Bx2g2+l"Months after initialization" -By0.2g0.2 >> $ps
for i in {1..10}
do
   gmt psxy  -R -J -N -W0.75p,$color1b -O -K rmse_tauv_${i}.txt >> $ps
done
gmt psxy  -R -J -N -W1.5p,$color1a -O -K rmse_tauv_mean.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2a -O -K rmse_tauv_ecco_obs_seas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K rmse_tauv_seas.txt >> $ps
echo "1 1.4 d" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O >> $ps


gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
