set term postscript color
set out "output.ps"
plot "0.05.txt" title "step=0.05" w l,"0.10.txt" title "step=0.10"  w l, "0.25.txt" title "step=0.25" w l, "0.50.txt" w l, "original.txt" w l
set title "Plots"
