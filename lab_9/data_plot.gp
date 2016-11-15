set term postscript color
set out "Heuns_method.ps"
set title "Euler-Heun's Comparison"
plot "0.50.txt" u 1:2 title "Exact" w lp, "0.50.txt" u 1:3 title "Heun's" w lp, "0.50.txt" u 1:5 title "Euler" w lp  

