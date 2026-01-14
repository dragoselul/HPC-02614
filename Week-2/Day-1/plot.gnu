#!/usr/bin/gnuplot
reset

set terminal pngcairo size 1200,900 enhanced font 'Arial,12'

# -------- Cache sizes (BYTES) - edit for YOUR CPU --------
L1 = 32.0*1024.0           # 32 KiB L1 data cache (typical)
L2 = 256.0*1024.0          # 256 KiB L2 cache (typical)
L3 = 8.0*1024.0*1024.0     # 8 MiB L3 cache (typical)
RAM = 16.0*1024.0*1024.0   # 16 MiB - beyond this = main RAM only

# Convert to KiB (data is already in KiB)
KiB = 1024.0
L1k = L1/KiB
L2k = L2/KiB
L3k = L3/KiB
RAMk = RAM/KiB

# -------- Shared plot styling --------
set xlabel "Working-set size (KiB)"
set ylabel "MFLOPS"
set logscale x 2
set grid
set key right bottom

# Cache boundary markers (vertical lines spanning full plot height)
set style line 101 lw 2 dt 2 lc rgb "#888888"
set style line 102 lw 2 dt 2 lc rgb "#666666"
set style line 103 lw 2 dt 2 lc rgb "#444444"
set style line 104 lw 2 dt 3 lc rgb "#cc0000"

set arrow 11 from L1k, graph 0 to L1k, graph 1 nohead ls 101
set arrow 12 from L2k, graph 0 to L2k, graph 1 nohead ls 102
set arrow 13 from L3k, graph 0 to L3k, graph 1 nohead ls 103
set arrow 14 from RAMk, graph 0 to RAMk, graph 1 nohead ls 104

set label 21 "L1" at L1k, graph 0.98 center font ",10"
set label 22 "L2" at L2k, graph 0.98 center font ",10"
set label 23 "L3" at L3k, graph 0.98 center font ",10"
set label 24 "RAM" at RAMk, graph 0.98 center font ",10"

# ============ Plot 1: All permutations ============
set output 'permutation_comparison.png'
set title "Matrix Multiplication Performance: Loop Order vs Memory (with cache boundaries)"

plot 'nat.dat' using 1:2 with linespoints lw 2 pt 7 title 'nat (ijk)', \
     'mnk.dat' using 1:2 with linespoints lw 2 pt 5 title 'mnk', \
     'mkn.dat' using 1:2 with linespoints lw 2 pt 9 title 'mkn', \
     'kmn.dat' using 1:2 with linespoints lw 2 pt 11 title 'kmn', \
     'knm.dat' using 1:2 with linespoints lw 2 pt 13 title 'knm', \
     'nmk.dat' using 1:2 with linespoints lw 2 pt 6 title 'nmk', \
     'nkm.dat' using 1:2 with linespoints lw 2 pt 8 title 'nkm'

# ============ Plot 2: Best vs BLAS ============
set output 'native_vs_blas.png'
set title "Native Implementation vs BLAS (with cache boundaries)"
set key left top

plot 'mkn.dat' using 1:2 with linespoints lw 2 pt 7 title 'mkn (best native)', \
     'lib.dat' using 1:2 with linespoints lw 2 pt 5 title 'BLAS (lib)'

# ============ Plot 3: Block size comparison ============
set output 'block_comparison.png'
set title "Blocked Matrix Multiplication: Block Size vs Memory (with cache boundaries)"
set key right bottom

plot 'blk_8.dat'   using 1:2 with linespoints lw 2 pt 7  title 'bs=8', \
     'blk_16.dat'  using 1:2 with linespoints lw 2 pt 5  title 'bs=16', \
     'blk_32.dat'  using 1:2 with linespoints lw 2 pt 9  title 'bs=32', \
     'blk_64.dat'  using 1:2 with linespoints lw 2 pt 11 title 'bs=64', \
     'blk_128.dat' using 1:2 with linespoints lw 2 pt 13 title 'bs=128', \
     'blk_256.dat' using 1:2 with linespoints lw 2 pt 6  title 'bs=256'

# ============ Plot 4: Multiplot ============
set output 'all_performance.png'
set terminal pngcairo size 1400,1000 enhanced font 'Arial,11'
set multiplot layout 2,2 title "Matrix Multiplication Performance (with cache + RAM boundaries)"

# Top left
set title "Loop Order Comparison"
set key right bottom font ",9"
plot 'nat.dat' using 1:2 with linespoints lw 2 pt 7 ps 0.7 title 'nat', \
     'mkn.dat' using 1:2 with linespoints lw 2 pt 5 ps 0.7 title 'mkn', \
     'kmn.dat' using 1:2 with linespoints lw 2 pt 9 ps 0.7 title 'kmn', \
     'knm.dat' using 1:2 with linespoints lw 2 pt 11 ps 0.7 title 'knm'

# Top right
set title "Native vs BLAS"
set key left top
plot 'mkn.dat' using 1:2 with linespoints lw 2 pt 7 title 'mkn', \
     'lib.dat' using 1:2 with linespoints lw 2 pt 5 title 'BLAS'

# Bottom left
set title "Block Size Effect"
set key right bottom
plot 'blk_16.dat' using 1:2 with linespoints lw 2 pt 7 title 'bs=16', \
     'blk_32.dat' using 1:2 with linespoints lw 2 pt 5 title 'bs=32', \
     'blk_64.dat' using 1:2 with linespoints lw 2 pt 9 title 'bs=64'

# Bottom right
set title "Best Blocked vs BLAS"
set key left top
plot 'blk_64.dat' using 1:2 with linespoints lw 2 pt 7 title 'blk (bs=64)', \
     'lib.dat'    using 1:2 with linespoints lw 2 pt 5 title 'BLAS'

unset multiplot

print "Plots generated with cache boundaries: L1, L2, L3, RAM"
