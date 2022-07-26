
rm *.ps
ps=corr_rmse_seas.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14.5cx9c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1.5p
gmt set MAP_TITLE_OFFSET=0.05c
gmt set MAP_TICK_LENGTH_PRIMARY=0.05c
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

Rrmse=1/12/4.4/9.4
Rcorr=1/12/-0.2/0.7
J=X3.3c/3.5c
Bx=x3g3
Byrmse=y1g1
Bycorr=y0.2g0.2

function print_rmse {
  gmt psxy  -R -J -N -W1.5p,$color5  -O -K rmse_ssh_$seas.txt >> $ps  
  gmt psxy  -R -J -N -W1.5p,$color4  -O -K rmse_pers_$seas.txt >> $ps  
  gmt psxy  -R -J -N -W1.5p,$color1a  -O -K rmse_pred_$seas.txt >> $ps
  gmt psxy  -R -J -N -W1.5p,$color2a  -O -K rmse_seas_$seas.txt >> $ps
  gmt psxy  -R -J -N -W1.5p,$color0  -O -K rmse_mean_$seas.txt >> $ps
}

function print_corr {
  gmt psxy  -R -J -N -W1.5p,$color5  -O -K corr_ssh_$seas.txt >> $ps
  gmt psxy  -R -J -N -W1.5p,$color4  -O -K corr_pers_$seas.txt >> $ps
  gmt psxy  -R -J -N -W1.5p,$color1a  -O -K corr_pred_$seas.txt >> $ps
  gmt psxy  -R -J -N -W1.5p,$color2a  -O -K corr_seas_$seas.txt >> $ps

}

seas="DJF"
gmt psbasemap -K -R$Rrmse -J$J -X0.7c -Y5.1c -BWesn+t"DJF" -B$Bx -B$Byrmse > $ps
print_rmse 
echo "1 6.3 RMSE (cm)" |  gmt pstext -D-0.6c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps

seas="MAM"
gmt psbasemap -O -K -R$Rrmse -J$J -X3.45c -Bwesn+t"MAM" -B$Bx -B$Byrmse >> $ps
print_rmse 

seas="JJA"
gmt psbasemap -O -K -R$Rrmse -J$J -X3.45c -Bwesn+t"JJA" -B$Bx -B$Byrmse >> $ps
print_rmse 

seas="SON"
gmt psbasemap -O -K -R$Rrmse -J$J -X3.45c -Bwesn+t"SON" -B$Bx -B$Byrmse >> $ps
print_rmse 

### Corr
seas="DJF"
gmt psbasemap -O -K -R$Rcorr -J$J -X-10.35c -Y-3.8c -BWeSn -B$Bx -B$Bycorr >> $ps
print_corr 
echo "1 0.25.4 ACC (-)" |  gmt pstext -D-0.6c/0c -R -J -N -F+f7,Hind-SemiBold+jCM+a90 -O -K >> $ps

seas="MAM"
gmt psbasemap -O -K -R$Rcorr -J$J -X3.45c -BweSn -B$Bx -B$Bycorr >> $ps
print_corr 

seas="JJA"
gmt psbasemap -O -K -R$Rcorr -J$J -X3.45c -BweSn -B$Bx -B$Bycorr >> $ps
print_corr 
echo "1 -0.2 Months after initialization" |  gmt pstext -D0.0c/-0.45c -R -J -N -F+f7,Hind-Light+jCM -O -K >> $ps

seas="SON"
gmt psbasemap -O -K -R$Rcorr -J$J -X3.45c -BweSn -B$Bx -B$Bycorr >> $ps
print_corr 

gmt psbasemap -O -K -R0/4.9/-0.5/0.5 -JX14.15c/0.5c -X-10.85c -Y-1.25c -T100/100/1 >> $ps


echo -e "0.0 0 \n 0.1 0" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.105 0 Mean seasonal cycle" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "1.04 0 \n 1.14 0" | gmt psxy -R -J -O -K -W1.5p,$color1a >> $ps
echo "1.145 0 Hybrid (ECCO-CCSM4)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps


echo -e "2.19 0 \n 2.29 0" | gmt psxy -R -J -O -K -W1.5p,$color2a >> $ps
echo "2.295 0 Hybrid (ECCO-DynP)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "3.23 0 \n 3.33 0" | gmt psxy -R -J -O -K -W1.5p,$color4 >> $ps
echo "3.335 0 Damped persistence" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo -e "4.28 0 \n 4.38 0" | gmt psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "4.385 0 CCSM4 SLA" | gmt pstext -R -J -F+f8+jLM -O  >> $ps



gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

