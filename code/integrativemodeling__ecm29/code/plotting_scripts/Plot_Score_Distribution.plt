reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10

set border lw 5 lc rgb "#484848"

set boxwidth 0.1 absolute
set style fill empty border rgb "#EB7262"
set style fill empty border rgb "#000080"

unset key 

set encoding iso_8859_1
set ylabel "Number of Models ({/Symbol \264}10^3)" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"
set xlabel "Score"     	      tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 57"

set yr [0:3] noreverse nowriteback
set xr [1535:1640] noreverse nowriteback

set pointsize 2

set xtics 1535,50,1635 border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"
#set format x "%.2f"

set ytics 0, 1,3 nomirror tc rgb "#484848"

set label "U-test p-value: 0.00\nCliff's {/:Italic d}: 0.13"  at graph 0.50, 0.95 right font 'Arial-Bold, 40' front tc rgb "#484848"

set key tc rgb "#484848"
set output "ECM29_Score_Distribution.pdf"
plot      "TS_HA.txt" usi 1:($2/1000) w histeps lw 10 lc rgb "#EB7262" title "Sample 1", \
	  "TS_HB.txt" usi 1:($2/1000) w histeps lw 10 lc rgb "#000080" title "Sample 2"


set output