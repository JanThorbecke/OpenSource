set terminal png transparent nocrop enhanced  size 500,420 
set output 'binary.1.png'

#set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
#set hidden3d offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set style data lines
#set ticslevel 0
#set title "Hidden line removal of explicit binary surfaces" 
#set xrange [ 1.00000 : 100.00000 ] noreverse nowriteback
#set yrange [ 1.00000 : 100.00000 ] noreverse nowriteback
set zrange [ 1.00000 : 200000000.00000 ] 
splot "corr2d.bin" binary array=100x100 format="%float32" endian=little

