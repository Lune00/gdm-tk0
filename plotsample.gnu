set terminal epslatex color colortext solid size 12.5cm,7.5cm standalone background rgb "#00003B"
set size ratio 1/1.61
set key spacing 2.5 width 1.8
set output 'plotsample.tex'
set size ratio -1

unset ylabel
unset xlabel
unset border
unset xtics
unset ytics
set yr [-0.1:1.2]

a=0.08
v = 1.8
set cbrange [-1:1]

sample = "sample00000.txt"
walls ="wall00000.txt"
bands = "bands00000.txt"

set palette rgb 33,13,10
#sample u 1:2:($4*a/v):($5*a/v):(sqrt( ($4/v)**2 + ($5/v)**2 ) *($4/abs($4)) ) w vectors lc palette lw 1 notitle,\

plot sample u 1:2:3 w circles lc rgb "#4CB3FF" fs transparent solid 1. noborder notitle,\
     bands u 1:2:($3*0.8) w circles lc rgb "#00003B" fs transparent solid 1. notitle,\
     walls u 1:2:3 w circles lw 0.1 lc rgb 'white' fs transparent solid 1. border lc rgb 'black' notitle,\
     sample u 1:2:($4*a/v):($5*a/v)w vectors lc '#f4f4f4' lw 1.2 notitle
