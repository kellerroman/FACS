set yrange [-0.1:1.1]
set y2tics

#plot "start.csv" u 1:3 w l \
#   , "sol.csv" u 1:3 w l \
#   , "output"  u 2:3 w l
#   , "sol.csv" u 1:2 w l axis x1y2 \

plot "sol.csv" u 1:3 w l \
   , "exact.csv"  u 2:3 w l \
   , "sol.csv" u 1:4 w l axis x1y2 \
   , "exact.csv"  u 2:5 w l axis x1y2
