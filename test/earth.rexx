/* macro for Planet */
view 1600 1200
'map(|ppmforge -night -seed 203 -stars 200 -saturation 250 -xsize 1600 -ysize 1200)'
image
generate
'map(|jpegtopnm map/earth-map.jpg|pnmgamma 1.5)'
light 404040 'z-20' y15 x0
planet 200 800 2 x35 z0 y110
generate
'map(|pnmscale -width 2048 -height 1024 map/white.pgm,full)'
'alpha(|jpegtopnm map/cloud-map.jpg|ppmtopgm)'
planet 202 800 0 x35 z0 y110
generate
'save(|pnmtojpeg -quality 90 > earth.jpg)'
exit
