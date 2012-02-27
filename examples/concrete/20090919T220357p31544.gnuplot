#!/usr/bin/env gnuplot
#
# created Sat Sep 19 22:06:39 2009 (20090919_22:06)
#
set term wxt 0 persist
set xlabel 'eps'
set grid
set datafile missing 'nan'
set ylabel 'sigma'
plot  "< bzcat 20090919T220357p31544.data.bz2" using 1:3 title '← sigma(eps)' with lines
