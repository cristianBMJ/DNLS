set output "density.png"
set xlabel "Cn^2"
set ylabel "Tiempo"
set autoscale xfix
set autoscale yfix
set autoscale cbfix
set cbrange [0:1]
plot 'DLS23.txt' matrix with image 
pause mouse
reset

