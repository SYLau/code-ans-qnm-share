set term wxt size 680, 480 enhanced font "Century,12"

set encoding utf8

set xtics font ",14"
set ytics font ",14"

set xlabel "𝑓 (Hz)" font "Cambria,14"
set ylabel "|𝐴_{in}|" font "Cambria,14"

set logscale y

plot \
"ingoing_amp.txt" using 1:2 notitle w l lc 8