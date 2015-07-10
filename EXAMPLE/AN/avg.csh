#!/bin/csh

awk '{print $1}' CNBSE/absspct > CNBSE/abs.energy
rm -f CNBSE/abs.tmp
touch CNBSE/abs.tmp

foreach i ( 1  2  )

foreach j ( 1 2 3 )

awk '{print $3}' CNBSE/absspct_N_.${i}_1s_${j} > CNBSE/abs.${i}.${j}

end
paste CNBSE/abs.energy CNBSE/abs.${i}.1 CNBSE/abs.${i}.2 CNBSE/abs.${i}.3 | awk '{print $1, ($2+$3+$4)/3}' > CNBSE/abs.avg.${i}

paste CNBSE/abs.tmp CNBSE/abs.avg.${i} > CNBSE/abs.tmp2
mv CNBSE/abs.tmp2 CNBSE/abs.tmp

end
# no division by 8 so we match the 1 unit cell zero-temp run
awk '{print $1, ($2+$4+$6+$8)}' CNBSE/abs.tmp > CNBSE/abs.avg

