#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as g.sh
# chmod +x g.sh
# ./g.sh
# checked in https://www.shellcheck.net/

 
# for all pgm files in this directory
#!/bin/bash
for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename "$file" .pgm)
  # convert  using ImageMagic
  convert "${b}".pgm -resize 600x600 "${b}".png
  echo "$file"
done

 
echo OK
# end

