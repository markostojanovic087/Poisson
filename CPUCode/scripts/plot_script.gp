set dgrid3d 30,30
set hidden3d
splot "out.dat" using 1:2:3 "%lf%lf%lf" with lines
