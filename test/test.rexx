/* macro for Planet */
reset 1000
view 800 600
'map(|pngtopnm map/test-map.png,full,40)'
light 202020 'z-20' y45 x35
planet 200 800 6 'z-20' y45 x15
generate
'save(|pnmtojpeg -quality 90 > test.jpg)'
exit
