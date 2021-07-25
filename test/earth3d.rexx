/* macro for Planet */
view 1280 1440 6
'map(|ppmforge -night -seed 203 -stars 200 -saturation 250 -xsize 1600 -ysize 1200)'
image
generate '-10'
'map(|jpegtopnm map/earth-map.jpg|pnmgamma 1.5)'
light 404040 'z-20' y15 x0
planet 200 800 2 x35 z0 y110
generate 6
'map(|pnmscale -width 2048 -height 1024 map/white.pgm,full)'
'alpha(|jpegtopnm map/cloud-map.jpg|ppmtopgm)'
planet 205 800 0 x35 z0 y110
generate 6
'save(|pnmtojpeg -quality 90 > earth3d.jpg)'
exit
