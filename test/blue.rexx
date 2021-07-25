/* macro for Planet */
view 1600 1200
'map(|ppmforge -night -stars 200 -saturation 250 -xsize 1600 -ysize 1200)'
image
generate
'map(|jpegtopnm map/blue-map.jpg,full,40)'
light 202020 'z-25' y45 x35
planet 140 800 4 'z-20' y45 x15
generate
'map(|pngtopnm map/brown-ring.png)'
'alpha(|pngtopnm -alpha map/brown-ring.png)'
ring 220 340 60
generate
'save(|pnmtojpeg -quality 90 > blue.jpg)'
exit
