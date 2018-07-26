#!/bin/csh

awk '{print $1}' CNBSE/xesspct_N_.0001_1s_01 > CNBSE/xes.energy
rm -f CNBSE/xes.tmp
touch CNBSE/xes.tmp

foreach i ( 1  2  )

foreach j ( 1 2 3 )

awk '{print $3}' CNBSE/xesspct_N_.000${i}_1s_0${j} > CNBSE/xes.${i}.${j}

end
paste CNBSE/xes.energy CNBSE/xes.${i}.1 CNBSE/xes.${i}.2 CNBSE/xes.${i}.3 | awk '{print $1, ($2+$3+$4)/3}' > CNBSE/xes.avg.${i}

paste CNBSE/xes.tmp CNBSE/xes.avg.${i} > CNBSE/xes.tmp2
mv CNBSE/xes.tmp2 CNBSE/xes.tmp

end
# no division by 8 so we match the 1 unit cell zero-temp run
awk '{print $1, ($2+$4+$6+$8)}' CNBSE/xes.tmp > CNBSE/xes.avg

