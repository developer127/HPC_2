set terminal svg size 900, 500
set output "runtime.svg"

# logarithmic scale for the y axes makes sense for exponentally growing data.
set logscale x
set logscale y

set xlabel "Matrix dim A: M=N" font ",16"
set ylabel "Time [ns]" font ",16"
set title "Solving linear system Ax=b" font ",18"
set key outside
set pointsize 0.5
plot "runtime.dat" every 14::8 using "%*30[^\n]%lf%*2[\n]%*30[^\n]%lf" with linespoints lt 2 lw 3 title "fulMatrix"
#, \
#    "benchmarks/bench_blas.data"\
#       using 1:4 with linespoints lt 3 lw 3 title "tripletFormat"
