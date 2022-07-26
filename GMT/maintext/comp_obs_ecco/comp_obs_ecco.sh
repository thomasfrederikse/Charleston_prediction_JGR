rm *.ps
ps=comp_obs_ecco.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_14.5cx9c
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

color0=40/40/40
color1=#4c78a8
color2=#f58518


gmt psbasemap -K -R1995/2017/-21/19 -JX13.2c/4.0c -X0.8c -Y4.7c -BWesn -Bx3g3 -By5g5+l"Sea level (cm)" > $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K alt.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2 -O -K ECCO_ssh.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1 -O -K ECCO_adj.txt >> $ps
echo -e "2009.5 -13 \n 2010.5 -13" | gmt psxy -R -J -O -K -W2p,$color0 >> $ps
echo "2010.6 -13 Observed SLA (Altimetry)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "2009.5 -16 \n 2010.5 -16" | gmt psxy -R -J -O -K -W2p,$color2 >> $ps
echo "2010.6 -16 (@~r@~ = 0.66) ECCO SLA" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "2009.5 -19 \n 2010.5 -19" | gmt psxy -R -J -O -K -W2p,$color1 >> $ps
echo "2010.6 -19 (@~r@~ = 0.66) ECCO reconstruction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo "1995 -21 a" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O -K >> $ps


gmt psbasemap -O -K -R1995/2017/-21/19 -JX13.2c/4.0c  -Y-4.4c -BWeSn -Bx3g3+l"Time (years)" -By5g5+l"Sea level (cm)" >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K tg.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color2 -O -K ECCO_ssh.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1 -O -K ECCO_adj.txt >> $ps
echo -e "2009.5 -13 \n 2010.5 -13" | gmt psxy -R -J -O -K -W2p,$color0 >> $ps
echo "2010.6 -13 Observed SLA (tide gauge)" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "2009.5 -16 \n 2010.5 -16" | gmt psxy -R -J -O -K -W2p,$color2 >> $ps
echo "2010.6 -16 (@~r@~ = 0.67) ECCO SLA" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "2009.5 -19 \n 2010.5 -19" | gmt psxy -R -J -O -K -W2p,$color1 >> $ps
echo "2010.6 -19 (@~r@~ = 0.64) ECCO reconstruction" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo "1995 -21 b" | gmt pstext -D0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLB -O  >> $ps

gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

