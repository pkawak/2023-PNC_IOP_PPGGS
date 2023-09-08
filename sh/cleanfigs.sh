#!/bin/bash

name=`echo $1 | sed 's/\.tex//'`

mkdir -p "unused"
for image_file in $(ls figs/*.*); do

    imgcount=`grep $image_file $name.log -c`
    if [[ "$imgcount" == "0" ]]; then
      echo "File $image_file is not in use."
      mv "$image_file" "unused/"
    else
      echo "File $image_file is in use."
    fi

done

