#!/bin/sh

#ai2nbse build/git directory
AI="/users/kgilmore/Code/AI2NBSE/"

# ocean install directory
OCEAN_INSTALL="/users/kgilmore/Code/OCEAN-develop/BIN/"

echo ${AI}

for file in bridgegw.x bridgelad.x conbbtoug.x conbbtoux.x condens2.x ll.x spect.x whom0.x
do

cp ${AI}/Src/NBSE/zbridge/zexe/$file $OCEAN_INSTALL/

done

for file in psiln.pl stepper.h
do
cp ${AI}/Src/NBSE/$file $OCEAN_INSTALL/

done
cp ${AI}/Src/NBSE/zexe/smnanl.x $OCEAN_INSTALL/

cp RIXS.pl $OCEAN_INSTALL/
cp rixs.h $OCEAN_INSTALL/
