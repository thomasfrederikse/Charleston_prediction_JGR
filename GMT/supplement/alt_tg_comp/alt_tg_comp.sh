rm *.ps
ps=alt_tg_comp.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_11.2cx10.9c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1.5p
gmt set MAP_TITLE_OFFSET=0.05c
gmt set MAP_TICK_LENGTH_PRIMARY=0.03c
gmt set MAP_TICK_PEN=thinnest,40,40,40
gmt set MAP_GRID_PEN_PRIMARY=thinnest,lightgrey
gmt set MAP_ANNOT_OFFSET_PRIMARY=2p
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,black
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,black
gmt set FONT_LABEL                     = 8p,Hind-Light,black
gmt set FONT_LOGO                      = 8p,Hind-Light,black
gmt set FONT_TITLE                     = 8p,Hind-Light,black

color0=60/60/60 # MSS
color0b=120/120/120 # MSS
color1a=#4c78a8 # ECCO rec
color1b=#9ecae9
color2a=#f58518 # Hybrid rec ens
color2b=#ffbf79 # Hybrid rec mean


gmt psbasemap -K -R1995/2017/-25/32 -JX5c/4c -X0.82c -Y6.5c -BWeSn+t"Seasonal cycle retained" -Bx5g5 -By8g8+l"Height (cm)" > $ps
gmt psxy  -R -J -N -W1.5p,$color0  -O -K tg.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1a -O -K alt.txt >> $ps
echo "2017 -25 a" | gmt pstext -D-0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jRB -O -K >> $ps

gmt psbasemap -O -K -R1995/2017/-25/32 -JX5c/4c -X5.25c -BweSn+t"Seasonal cycle removed" -Bx5g5 -By8g8 >> $ps
gmt psxy  -R -J -N -W1.5p,$color0  -O -K tg_noseas.txt >> $ps
gmt psxy  -R -J -N -W1.5p,$color1a -O -K  alt_noseas.txt >> $ps
echo -e "2008 29 \n 2010.5 29" | gmt psxy -R -J -O -K -W1.5p,$color1a >> $ps
echo "2010.9 29 Altimetry" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps
echo -e "2008 24 \n 2010.5 24" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "2010.9 24 Tide gauge" | gmt pstext -R -J -F+f8+jLM -O -K >> $ps

echo "2017 -25 b" | gmt pstext -D-0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jRB -O -K >> $ps

gmt psbasemap -O -K -R-25/30/-25/30 -JX5c/5c -X-5.25c -Y-5.8c -BWeSn -Bx8g8+l"Altimetry (cm)" -By8g8+l"Tide gauge (cm)" >> $ps
gmt psxy  -R -J -N -W1.5p,$color0b  -O -K straightline.txt >> $ps
gmt psxy  -R -J -N -O -K -Sc0.15c -G$color1a tg_vs_alt.txt >> $ps
echo "30 -25 c" | gmt pstext -D-0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jRB -O -K >> $ps

gmt psbasemap -O -K -R-25/30/-25/30 -JX5c/5c -X5.25c -BweSn -Bx8g8+l"Altimetry (cm)" -By8g8 >> $ps
gmt psxy  -R -J -N -W1.5p,$color0b  -O -K straightline.txt >> $ps
gmt psxy  -R -J -N -O -K -Sc0.15c -G$color1a tg_vs_alt_noseas.txt >> $ps
echo "30 -25 d" | gmt pstext -D-0.125c/0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jRB -O >> $ps

gmt psconvert -Tg -E400 -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
