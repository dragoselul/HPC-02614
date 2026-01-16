# plot_scaling.gnu - Plot thread scaling performance
#
# Usage: gnuplot plot_scaling.gnu
#

set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'thread_scaling.png'

set multiplot layout 2,2 title "OpenMP Thread Scaling Analysis" font ',14'

# ============================================
# Plot 1: GFLOPS vs Threads (Actual vs Ideal)
# ============================================
set title "Performance Scaling"
set xlabel "Number of Threads"
set ylabel "GFLOPS"
set grid
set key top left

plot 'thread_scaling.dat' using 1:3 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Actual Performance', \
     'thread_scaling.dat' using 1:4 with lines lw 2 lc rgb '#dd181f' dashtype 2 title 'Ideal (Linear) Scaling'

# ============================================
# Plot 2: Speedup vs Threads
# ============================================
set title "Speedup"
set xlabel "Number of Threads"
set ylabel "Speedup (x)"
set key top left

# Get max threads from data for ideal line
stats 'thread_scaling.dat' using 1 nooutput
max_threads = STATS_max

plot 'thread_scaling.dat' using 1:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Actual Speedup', \
     'thread_scaling.dat' using 1:1 with lines lw 2 lc rgb '#dd181f' dashtype 2 title 'Ideal Speedup'

# ============================================
# Plot 3: Efficiency vs Threads
# ============================================
set title "Parallel Efficiency"
set xlabel "Number of Threads"
set ylabel "Efficiency (%)"
set yrange [0:110]
set key top right

plot 'thread_scaling.dat' using 1:6 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Efficiency', \
     100 with lines lw 2 lc rgb '#dd181f' dashtype 2 title 'Ideal (100%)'

# ============================================
# Plot 4: Execution Time vs Threads
# ============================================
set title "Execution Time"
set xlabel "Number of Threads"
set ylabel "Time (seconds)"
set yrange [*:*]
set key top right

# Get baseline time (single thread) from data
stats 'thread_scaling.dat' using 2 every ::0::0 nooutput
t1 = STATS_min

plot 'thread_scaling.dat' using 1:2 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Actual Time', \
     'thread_scaling.dat' using 1:(t1/$1) with lines lw 2 lc rgb '#dd181f' dashtype 2 title 'Ideal Time'

unset multiplot

# Reset terminal for additional single plots
set terminal pngcairo size 800,600 enhanced font 'Arial,12'

# ============================================
# Single plot: Performance comparison
# ============================================
set output 'performance_scaling.png'
set title "Jacobi Solver: Actual vs Ideal Performance" font ',14'
set xlabel "Number of Threads"
set ylabel "GFLOPS"
set grid
set key top left

plot 'thread_scaling.dat' using 1:3 with linespoints pt 7 ps 2 lw 2 lc rgb '#0060ad' title 'Actual Performance', \
     'thread_scaling.dat' using 1:4 with lines lw 2 lc rgb '#dd181f' dashtype 2 title 'Ideal (Linear) Scaling'

print "Plots generated: thread_scaling.png, performance_scaling.png"

