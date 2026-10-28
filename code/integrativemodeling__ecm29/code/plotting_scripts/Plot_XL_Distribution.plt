reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10

set border lw 5 lc rgb "#484848"

set boxwidth 0.1 absolute
set style fill empty border rgb "#EB7262"
set style fill empty border rgb "#8666FB"

unset key 

set encoding iso_8859_1
set ylabel "Number of Cross-links (AU)" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"
set xlabel "Euclidian distance ({/Symbol 305})" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"

set yr [0:1] noreverse nowriteback
set xr [0:100] noreverse nowriteback

set pointsize 2

set xtics 0,25,100 border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"
set format y "%.2f"

set ytics 0,.25,1 nomirror tc rgb "#484848"

set linetype 5 dashtype 2 lw 10
set arrow nohead from 30,0 to 30,1 lc rgb "red" lt 5

set key tc rgb "#484848"
set output "ECM29_XL_C1.pdf"
plot      "XLs_H.txt" usi 1:($2/50000) w histeps lw 10 lc rgb "#484848" notitle


set output