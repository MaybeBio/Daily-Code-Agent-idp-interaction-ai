names = "ECM29"

do for [file in names]{ 
reset

set terminal pdfcairo enhanced color font "Arial-Bold,40" size 10,10
set border lw 3 lc rgb "#484848"

set boxwidth 0.5 absolute
set style fill solid 1.00 border rgb "#484848"
unset key 

set encoding iso_8859_1
unset xlabel

set ylabel  "Number of models ({/Symbol \264}10^3)" tc rgb "#484848" offset 0,0 font "Arial-Bold, 57"
set y2label "Population (%)" tc rgb "#484848" offset 0,0 font "Arial-Bold, 57"
set format y "%.1f"

set pointsize 2
set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set xtics  norangelimit
set xtics   ()
set ytics 0,2.5,7.5 nomirror
set y2tics 0,20,100 nomirror
set label "{/Symbol c}^2-test p-value: 0.10\nCramer's V: 0.01" at graph 0.95, 0.96 right font 'Arial-Bold, 40' front tc rgb "#484848"
set xr [0:3.25]
set yr [0:7.5]
set y2r [0:75]

set output sprintf("%s_Population.pdf", file)
plot sprintf("%s.P.txt", file) usi ($2-0.5):($3/1000):xtic(sprintf("Cluster %i", $1)) w boxes notitle lc rgb "#EB7262", \
     "" 		       usi 2:($4/1000) w boxes notitle lc "#000080"

set output
}

