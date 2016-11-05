#! /usr/bin/env bash

export LC_ALL=C

echo "Combaining files..."

mark="total_"
SCALE="1.0"

#special treatment
#if [ -e mcdif_06082016 ]; then SCALE="200.0"; fi
#tmp_pth=`echo $PWD | grep -o 'S[0-9]' | tail -1`; 
#tmp_pth=`wc -l < ../jobs_list`
#if [ $tmp_pth -eq 1 ]; then SCALE="0.2"; fi
#echo $tmp_pth $SCALE


#last_step=0
#last_time=0
#mask=N.dat
#rm ${mark}${mask}
#for i in ./*[0-9]$mask; do
#    awk -v step=$last_step -v time=$last_time '{printf "%f %f", ($1+step),($2+time); for (i=3;i<=NF;i++){printf " %s", $i;};printf "\n";}' $i >> ${mark}${mask}
#    last_step=`awk -v st=$last_step 'END{print $1+st}' $i`
#    last_time=`awk -v st=$last_time 'END{print $2+st}' $i`
#done

#last_step=0
#last_time=0
#mask=E.dat
#rm ${mark}${mask}
#for i in ./*[0-9]$mask; do
#    awk -v step=$last_step -v time=$last_time '{printf "%f %f", ($1+step),($2+time); for (i=3;i<=NF;i++){printf " %s", $i;};printf "\n";}' $i >> ${mark}${mask}
#    last_step=`awk -v st=$last_step 'END{print $1+st}' $i`
#    last_time=`awk -v st=$last_time 'END{print $2+st}' $i`
#done

#if [-e error_log] 

last_step=0
last_time=0
mask=dR2.dat
rm -f ${mark}${mask}
rm -f files 
for i in *[0-9]$mask; do
    echo ${i%%${mask}} >> files
done
sort -n files > tmp
mv tmp files


for i in $(cat files); do
    FILE=$i${mask}
    if ! [ -e bad_data ]; then
    awk -v k=$SCALE -v step=$last_step -v time=$last_time '{if($2 == "inf"){print "bad data" > "bad_data";exit;} printf "%f %f", ($1+step),($2*k+time); for (i=3;i<=NF;i++){printf " %s", $i;};printf "\n";}' $FILE >> ${mark}${mask}
    last_step=`awk 'END{print $1}' ${mark}${mask}`
    last_time=`awk 'END{print $2}' ${mark}${mask}`
    fi
done

last_step=0
last_time=0
mask=dR.dat
rm -f ${mark}${mask}
for i in $(cat files); do
    FILE=$i${mask}
    if ! [ -e bad_data ]; then
    awk -v k=$SCALE -v step=$last_step -v time=$last_time '{if($2 == "inf"){print "bad data" > "bad_data";exit;} printf "%f %f", ($1+step),($2*k+time); for (i=3;i<=NF;i++){printf " %s", $i;};printf "\n";}' $FILE >> ${mark}${mask}
    last_step=`awk 'END{print $1}' ${mark}${mask}`
    last_time=`awk 'END{print $2}' ${mark}${mask}`
    fi
done
datar2="${mark}dR2.dat"
datar="${mark}dR.dat"

if ! [ -e bad_data ]; then
sort -n -k2,2 $datar2 > tmp && mv tmp $datar2
sort -n -k2,2 $datar > tmp && mv tmp $datar

awk '{print $7, $8, $13, $14, $19, $20, $25, $26}' $datar > r1.tmp
awk '{print $6, $10, $14, $18}' $datar2 > r2.tmp
#paste r1.tmp r2.tmp > r.tmp
#awk '{print $1,$4,$2,$5,$3,$6}' r.tmp > jumps.tmp
paste -d " " $datar2 r1.tmp > file1.tmp
mv file1.tmp ${datar2}
paste -d " " $datar r2.tmp > file2.tmp
mv file2.tmp ${datar}

rm *.tmp

#delete last line
head -n -1 $datar > tmp; mv tmp $datar;
head -n -1 $datar2 > tmp; mv tmp $datar2;
fi
#export notification
if [ -e bad_data ]; then
    rm ${mark}*
    echo $PWD >> $1
fi

