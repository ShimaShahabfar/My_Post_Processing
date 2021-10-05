#!/bin/bash 

if [ -f TMP0 ]
then
    rm TMP0
fi

if [ -f TMP1 ]
then
    rm TMP1
fi

if [ -f TMP2 ]
then
    rm TMP2
fi

if [ -f TMP3 ]
then
    rm TMP3
fi

if [ -f BaderCharge_Analysis.dat ]
then
    rm BaderCharge_Analysis.dat
fi
  

printf "%5s %5s %6s\n"  "#  Num" "  Atom" "  Bader"  > BaderCharge_Analysis.dat
echo "-------------------------" >> BaderCharge_Analysis.dat

grep -v '#\|\-\|VACUUM\|VACUUM\|NUMBER' ACF.dat | awk '{print $5}' > TMP1


ncol=$(grep -B 2 -i 'direct' POSCAR  | head -n 1  | awk '{print NF}')

count=0
for i in $(seq 1 1 $ncol);
do
    atm=$(grep -B 2 -i 'direct' POSCAR  | head -n 1  | awk -v var=$i '{print $var}')
    num=$(grep -B 2 -i 'direct' POSCAR  | tail -n 2  | head -n 1  | awk -v var=$i '{print $var}')
    chg=$(grep -A 1 "PAW_PBE $atm" POTCAR | grep -v 'TITEL\|LULTRA' | tail -n 2 | head -n 1)
    for j in $(seq 1 1 $num)
    do
        count=$(($count+1))
        printf "%5s %5s %6s\n"  "$count" "$atm" "$chg"  >> TMP0
    done
done

paste TMP0  TMP1 > TMP2


awk '{printf "%5d  %5s %7.2f\n", $1, $2, $3-$4}' TMP2 >> BaderCharge_Analysis.dat

grep -v -w "#\|\-\-" BaderCharge_Analysis.dat | sort -u -k3,3 | sort -k1,1 -n  > TMP3
echo "" >>  BaderCharge_Analysis.dat 
echo "------UNIQUE CHARGES-------" >>  BaderCharge_Analysis.dat 
echo "" >>  BaderCharge_Analysis.dat 

cat TMP3 >> BaderCharge_Analysis.dat
echo "----------------------------" >>  BaderCharge_Analysis.dat 

rm TMP0 TMP1 TMP2 TMP3
