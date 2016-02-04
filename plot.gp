set size square
set xrange [0:L-1]
set yrange [0:L-1]
set cbrange[-1:1]
unset key

load '../colors.gp'

set palette negative
do for [i=1000:250000:200] {
	set title sprintf("t = %d",i)
	plot './data/aout_'.i matrix w image, './data/tout_'.i ps 4 lc 2
	pause 0.03
}
