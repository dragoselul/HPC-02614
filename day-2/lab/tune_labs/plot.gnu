#!/usr/bin/gnuplot

set terminal pngcairo size 1200,1000 enhanced font 'Arial,12'
set output 'performance_comparison.png'

set multiplot layout 2,2 title "AOS vs SOA Performance Comparison"

# Plot 1: Mflop/s distance
set title "Performance: distance()"
set xlabel "Memory Footprint (kB)"
set ylabel "Mflop/s"
set logscale x
set grid
set key left top
plot 'aos.gcc.dat' using 1:2 with linespoints lw 2 pt 7 title 'AOS', \
     'soa.gcc.dat' using 1:2 with linespoints lw 2 pt 5 title 'SOA'

# Plot 2: Mflop/s total
set title "Performance: Total"
set xlabel "Memory Footprint (kB)"
set ylabel "Mflop/s"
plot 'aos.gcc.dat' using 1:4 with linespoints lw 2 pt 7 title 'AOS', \
     'soa.gcc.dat' using 1:4 with linespoints lw 2 pt 5 title 'SOA'

# Plot 3: Runtime
set title "Runtime Comparison"
set xlabel "Memory Footprint (kB)"
set ylabel "Runtime (s)"
set logscale y
plot 'aos.gcc.dat' using 1:5 with linespoints lw 2 pt 7 title 'AOS', \
     'soa.gcc.dat' using 1:5 with linespoints lw 2 pt 5 title 'SOA'

# Plot 4: Speedup
set title "AOS vs SOA Speedup"
set xlabel "Memory Footprint (kB)"
set ylabel "Speedup (AOS/SOA)"
unset logscale y
set key right top
plot 'aos.gcc.dat' using 1:($4/(column(4))) with linespoints lw 2 pt 7 notitle, \
     'soa.gcc.dat' using 1:($4/(column(4))) every ::0::0 with lines dt 2 lc 'gray' title 'Equal performance'

# Create individual plots as well
unset multiplot

set output 'mflops_distance.png'
set title "Performance: distance()"
set ylabel "Mflop/s"
set logscale x
unset logscale y
set key left top
plot 'aos.gcc.dat' using 1:2 with linespoints lw 2 pt 7 title 'AOS', \
     'soa.gcc.dat' using 1:2 with linespoints lw 2 pt 5 title 'SOA'

set output 'runtime.png'
set title "Runtime Comparison"
set ylabel "Runtime (s)"
set logscale xy
plot 'aos.gcc.dat' using 1:5 with linespoints lw 2 pt 7 title 'AOS', \
     'soa.gcc.dat' using 1:5 with linespoints lw 2 pt 5 title 'SOA'

print "Plots saved: performance_comparison.png, mflops_distance.png, runtime.png"
