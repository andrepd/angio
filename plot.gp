set size square
set xrange [0:127]
set yrange [0:127]
set cbrange[-1:1]
unset key

load '/home/andrepd/Development/gnuplot-colorbrewer/diverging/RdYlBu.plt'

set palette negative
do for [i=1000:100000:500] {
	set title sprintf("t = %d",i)
	plot './data/aout_'.i matrix w image, './data/tout_'.i ps 4
	pause 0.1
}
