#!/bin/bash
 
# script file for BASH 
# which bash
# save this file as g.sh
# chmod +x g.sh
# ./g.sh
# checked in https://www.shellcheck.net/
# shellcheck ./g.sh



echo "make all pgm files"
gcc d.c -lm -Wall -march=native -fopenmp
time ./a.out > a.txt

echo "convert all pgm files to png using Image Magic convert"
# for all pgm files in this directory
for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename "$file" .pgm)
  # convert  using ImageMagic
  convert "${b}".pgm -resize 2000x2000 "${b}".png
  echo "$file"
done



echo  "delete all pgm files"
rm ./*.pgm


# display OpenMP info
export  OMP_DISPLAY_ENV="TRUE"	

echo OK
# end
