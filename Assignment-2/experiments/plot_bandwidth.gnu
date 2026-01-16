# plot_bandwidth.gnu - Plot updates per second and memory bandwidth
#
# Usage: gnuplot plot_bandwidth.gnu
#
# Memory bandwidth estimation for 3D Jacobi stencil:
#   Per point update: 7 reads (6 neighbors + f) + 1 write = 8 doubles = 64 bytes
#   With cache line effects, effective bandwidth may be higher (~72-96 bytes/update)

set terminal pngcairo size 1400,600 enhanced font 'Arial,12'
set output 'bandwidth_analysis.png'

set multiplot layout 1,2 title "Poisson 3D Solver: Memory Bandwidth Analysis" font ',14'

# Bytes per grid point update (7 reads + 1 write = 8 doubles)
BYTES_PER_UPDATE = 64.0

# Peak memory bandwidth of your system (adjust for your hardware)
# e.g., DDR4-2666: ~85 GB/s, DDR4-3200: ~100 GB/s
PEAK_BW_GBS = 100.0

# ============================================
# Plot 1: Updates per Second vs Grid Size
# ============================================
set title "Grid Point Updates per Second"
set xlabel "Grid Size (N)"
set ylabel "Million Updates/s"
set format y "%.0f"
set grid
set key top right

# Data is in updates/s, convert to millions for display
plot 'updates_per_second.dat' using 1:($2/1.0e6) with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'MUpdates/s'

# ============================================
# Plot 2: Memory Bandwidth vs Grid Size
# ============================================
set title "Effective Memory Bandwidth"
set xlabel "Grid Size (N)"
set ylabel "Bandwidth (GB/s)"
set format y "%.1f"
set grid
set key top right

# Bandwidth (GB/s) = Updates/s * bytes_per_update / 1e9
plot 'updates_per_second.dat' using 1:($2 * BYTES_PER_UPDATE / 1.0e9) with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Effective Bandwidth', \
     PEAK_BW_GBS with lines lw 2 lc rgb '#dd181f' dashtype 2 title sprintf('Peak (~%.0f GB/s)', PEAK_BW_GBS)

unset multiplot

# ============================================
# Detailed single plot: Bandwidth with cache regions
# ============================================
set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'memory_bandwidth.png'

# Cache sizes (adjust for your CPU)
L1_KB = 512
L2_KB = 8192
L3_MB = 32

# N values at cache boundaries
N_L1 = (L1_KB * 1024.0 / 24.0) ** (1.0/3.0) - 2
N_L2 = (L2_KB * 1024.0 / 24.0) ** (1.0/3.0) - 2
N_L3 = (L3_MB * 1024.0 * 1024.0 / 24.0) ** (1.0/3.0) - 2

set title "Memory Bandwidth vs Grid Size (Jacobi 3D Stencil)" font ',14'
set xlabel "Grid Size (N)"
set ylabel "Bandwidth (GB/s)"
set format y "%.1f"
set grid
set key top right

# Cache boundary lines
set arrow 1 from N_L1, graph 0 to N_L1, graph 1 nohead lc rgb '#ff0000' lw 2 dt 2
set arrow 2 from N_L2, graph 0 to N_L2, graph 1 nohead lc rgb '#00aa00' lw 2 dt 2
set arrow 3 from N_L3, graph 0 to N_L3, graph 1 nohead lc rgb '#0000ff' lw 2 dt 2

set label 1 "L1" at N_L1, graph 0.95 center tc rgb '#ff0000' font ',10'
set label 2 "L2" at N_L2, graph 0.95 center tc rgb '#00aa00' font ',10'
set label 3 "L3" at N_L3, graph 0.95 center tc rgb '#0000ff' font ',10'

plot 'updates_per_second.dat' using 1:($2 * BYTES_PER_UPDATE / 1.0e9) with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Effective Bandwidth', \
     PEAK_BW_GBS with lines lw 2 lc rgb '#dd181f' dashtype 2 title sprintf('Peak Memory BW (~%.0f GB/s)', PEAK_BW_GBS)

# ============================================
# Combined plot: GFLOPS and Bandwidth
# ============================================
set output 'performance_bandwidth.png'

set title "Performance and Memory Bandwidth vs Grid Size" font ',14'
set xlabel "Grid Size (N)"
set ylabel "GFLOPS" tc rgb '#0060ad'
set y2label "Bandwidth (GB/s)" tc rgb '#dd181f'
set y2tics
set ytics nomirror
set format y "%.1f"
set format y2 "%.1f"
set grid
set key top center

# Clear previous arrows/labels
unset arrow 1
unset arrow 2
unset arrow 3
unset label 1
unset label 2
unset label 3

# Redraw cache boundaries
set arrow 1 from N_L1, graph 0 to N_L1, graph 1 nohead lc rgb '#888888' lw 1 dt 3
set arrow 2 from N_L2, graph 0 to N_L2, graph 1 nohead lc rgb '#888888' lw 1 dt 3
set arrow 3 from N_L3, graph 0 to N_L3, graph 1 nohead lc rgb '#888888' lw 1 dt 3

set label 1 "L1" at N_L1, graph 0.05 center tc rgb '#888888' font ',9'
set label 2 "L2" at N_L2, graph 0.05 center tc rgb '#888888' font ',9'
set label 3 "L3" at N_L3, graph 0.05 center tc rgb '#888888' font ',9'

plot 'grid_scaling.dat' using 1:4 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'GFLOPS' axes x1y1, \
     'updates_per_second.dat' using 1:($2 * BYTES_PER_UPDATE / 1.0e9) with linespoints pt 5 ps 1.5 lw 2 lc rgb '#dd181f' title 'Bandwidth (GB/s)' axes x1y2

print "Plots generated: bandwidth_analysis.png, memory_bandwidth.png, performance_bandwidth.png"
