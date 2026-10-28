reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 15,10

set border lw 5 lc rgb "#545454"
set boxwidth 0.25 absolute

set style fill solid 1.0  border rgb "#545454"
#unset key 

set encoding iso_8859_1
set xlabel "Cross-link" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 48"
set ylabel "Distance ({\305})" tc rgb "#484848" offset 0.0,0.0 font "Arial-Bold, 48"

set yr [0:100] noreverse nowriteback
set xr [0:66] noreverse nowriteback

set pointsize 2

set xtics border in scale 0,0 mirror rotate 90 offset character 0, 0, 0 autojustify tc rgb "#484848" font "Arial-Bold, 12"
#set format y "%.2f"

set ytics 0,25,100 nomirror tc rgb "#484848"

set linetype 5 dashtype 2 lw 10

set arrow nohead from 1,35 to 56,35 lc rgb "#646464" lt 5
set arrow nohead from 59,35 to 66,35 lc rgb "#646464" lt 5

set arrow nohead from 57.0,0 to 57.0,100 lc rgb "#545454" lw 5
set arrow nohead from 58.0,0 to 58.0,100 lc rgb "#545454" lw 5

#set arrow nohead from 57,0   to 59,0 lc rgb "#FFFF00" lw 5 lt 1 
#set arrow nohead from 57,100 to 59,100 lc rgb "#FFFF00" lw 5 lt 1

set key at graph 0.50, 0.95 tc rgb "#545454"
set output "ECM29_XL_Stats.pdf"
plot  "XLs_Ids_Statistics_C0_C1_2.txt" usi 1:3:4:xtic(2) w errorbars pt 4 ps 1.5 lw 3 lc rgb "#1B9E77" title "Cluster 1",\
      ""			       usi 1:5:6:xtic(2) w errorbars pt 7 ps 1.3 lw 3 lc rgb "#D95D02" title "Cluster 2"


set output