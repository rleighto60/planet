/* macro for Planet */
view 1280 1440 4
'map(|ppmforge -night -seed 203 -stars 200 -saturation 250 -xsize 1280 -ysize 1440)'
image
generate '-20'
light 202020 'z-20' y45 x40
'map(|pngtopnm map/test-map.png,full,40)'
planet 200 800 6 'z-20' y45 x15
generate 6
'map(map/white.pgm)'
ring 300 450 100
generate 6
planet 200 800 0 'z-20' 'y-30' x15
moon 15 500 2
generate 6
'save(|pnmtojpeg -quality 90 > test.jpg)'
exit
