/* macro for Planet */
reset 1400
view 1600 1200
'map(|ppmforge -night -stars 200 -saturation 250 -xsize 1600 -ysize 1200)'
image
generate
'map(|jpegtopnm map/blue-map.jpg,full,40)'
light 202020 'z-25' y45 x35
planet 140 800 6 'z-20' y45 x15
generate
'map(|jpegtopnm map/brown-ring.jpg,full,0,100,010101)'
light 404040 'z-25' y45 x35
ring 220 340 100
generate
'save(|pnmtojpeg -quality 90 > blue.jpg)'
exit
