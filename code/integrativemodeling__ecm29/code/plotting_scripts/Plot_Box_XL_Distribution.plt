reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10

set border lw 5 lc rgb "#484848"
set boxwidth 5.0 absolute
set style fill solid 1.0  border rgb "#545454"

stats 'XLs_Integer_H2.txt' usi 2 name "A" nooutput

unset key 
set encoding iso_8859_1
set ylabel "Number of Cross-links" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"
set xlabel "Euclidian Distance ({\305})" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"

set yr [0:1.05] noreverse nowriteback
set xr [0:100] noreverse nowriteback

set pointsize 2

set xtics 0,25,100 border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"
set format y "%.2f"
set ytics 0,.25,1 nomirror tc rgb "#484848"

set linetype 5 dashtype 2 lw 10
set arrow nohead from 30,0 to 30,1.05 lc rgb "red" lt 5

set key tc rgb "#484848"
set output "ECM29_BOX_XL_C1.pdf"

plot  "XLs_Integer_H2.txt" usi 1:($2 / A_max) w boxes lc rgb "#545454" notitle

set output