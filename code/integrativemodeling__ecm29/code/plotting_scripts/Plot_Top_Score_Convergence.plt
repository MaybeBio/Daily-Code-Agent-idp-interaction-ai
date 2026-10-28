reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10
set border lw 5 lc rgb "#484848"

set boxwidth 0.1 absolute
set style fill empty border rgb "#EB7262"
set style fill empty border rgb "#8666FB"

unset key 

set encoding iso_8859_1
set xlabel "Number of Models ({/Symbol \264}10^4)" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"
set ylabel "Score"     	      tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"

set xr [0:11] noreverse nowriteback
set yr [1535:1640] noreverse nowriteback

set pointsize 2

set ytics 1535,50,1635 border mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"
set xtics 0,2.5,10 border mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"
set format x "%.1f"

set linetype 5 dashtype 2 lw 10
set arrow nohead from 0,1635 to 10,1635 lc rgb "red" lt 5

#set label "{/:Italic t}-test p-value: 0.00\nCohen's {/:Italic d}: 0.18"  at graph 0.50, 0.95 right font 'Arial-Bold, 40' front tc rgb "#484848"

set key tc rgb "#484848"
set output "ECM29_TopScoring.pdf"
plot "TopScoring.txt" usi ($1/10000):2:3 w errorbar pt 7 ps 1.5 lc rgb "#484848" notitle


set output