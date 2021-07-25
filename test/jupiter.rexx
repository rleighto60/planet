/* macro for Planet */
address planet
'view(800,600)'
'map(|ppmforge -night -stars 200 -saturation 250 -xsize 1600 -ysize 1200)'
'image'
'generate'
'map(|jpegtopnm map/jupiter-map.jpg)'
'light(202020,z-20,y45,x35)'
'planet(240,1000,2,z-20,x15,y20)'
'generate'
'save(|pnmtojpeg -quality 90 > jupiter.jpg)';
exit
