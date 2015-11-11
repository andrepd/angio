set term wxt enhanced

set pm3d map
set hidden3d
set xlabel "y"
set ylabel "x"

do for[i=1:20] {
	name=sprintf("aout_%d",i*5000)
	spl name matrix
	pause 1
}
