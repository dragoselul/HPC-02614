#!/bin/bash

gnuplot plot_gflops_vs_grid.gnu
gnuplot plot_bw_vs_grid.gnu
gnuplot plot_updates_vs_grid.gnu
#gnuplot plot_speedup_dual.gnu

echo "All plots generated."
