#!/bin/sh

 if [[ $# != 1 ]]; then
   echo usage: $0 your.structure.POSCAR.vasp
   exit
 fi

 filename=$1

 if [ ! -f $filename ]; then
   echo $filename does not exits.
   exit
 fi
 
 bohr=0.529177
 acell=`sed -n '2p' $filename|awk '{print $1}'`
 acell=`echo "scale=6;$acell / $bohr"|bc`

 echo acell { $acell $acell $acell }
 echo
 echo rprim {
 sed -n '3,5p' $filename
 echo }
 echo

 . ~/shirley_QE4.3/scripts/arvid/Bibliothek/periodic.table pbe

 elem=`sed -n '6p' $filename`

 znucl=""

 shopt -s nocasematch

 echo pp_list {
 for e in $elem; do
   z=`get_atomicnumber $e`
   echo $e.fhi
   if [[ "$z" = *[0-9]* ]]; then
     znucl="$znucl $z"
   else
     znucl="$znucl ?"
   fi
 done
 echo }
 echo

 shopt -u nocasematch

 ntypat=`echo $elem|awk '{print NF}'`
 echo ntypat $ntypat
 echo zsymb { $elem }
 echo znucl { $znucl }
 echo

 num=`sed -n '7p' $filename`

 natom=0
 typat=""
 for i in `seq 1 $ntypat`; do
   numi=`echo $num|awk -v i=$i '{print $i}'`
   for j in `seq 1 $numi`; do
     typat="$typat $i"
   done
   natom=$((natom+numi))
 done

 echo natom $natom
 echo typat { $typat }
 echo

 echo xred {
 sed -n "9,$((8+natom))p" $filename
 echo }
 
