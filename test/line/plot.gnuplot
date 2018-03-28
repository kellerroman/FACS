set yrange [-0.1:1.1]
set y2tics

plot "old.csv" u 1:3 w l \
   , "sol.csv" u 1:3 w l \
   , "sol.csv" u 1:2 w l axis x1y2
