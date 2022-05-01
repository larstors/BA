
set term eps
convertedData='./convertedData.cfg'
latticeData=sprintf("%d.hc", L);
outFile=sprintf("%d_hc.eps", L);

outFile=sprintf("%s.eps", file);

set output outFile
unset key  
set size ratio -cos(pi/6)
set xtics
set ytics
set tics format ""
set autos
#set xrange[0:L*cos(pi/6.0)]
set yrange[-4:(3*L/2.+1)/1.]
#set yrange[0:45]

unset colorbox

plot latticeData u 1:2 w l lc 0,\
     convertedData u 1:2:3 palette pt 7 ps 0.6 notit

        
 
