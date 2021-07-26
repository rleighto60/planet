/* macro for Planet */
view 1600 1200 4
'map(|pngtopnm map/test-map.png,full,40)'
light 202020 'z-20' y45 x40
planet 200 800 6 'z-20' y45 x15
generate 6
'map(map/white.pgm)'
ring 300 450 100
generate 6
planet 200 800 0 'z-20' 'y-30' x15
moon 15 500 6
generate 6
'save(|pnmtojpeg -quality 90 > test.jpg)'
exit
