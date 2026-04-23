#!/bin/bash

set -e

# ******************************************************************************
# This script generates a list of files used by ffmpeg to generate a movie.
# Assuming the list of files is saved in a file called "ffmpeg_list.txt", the
# following command can be used to generate a movie from it:
#
#    ffmpeg -f concat -safe 0 -i ffmpeg_list.txt -c:v libx264 -crf 17 -preset veryslow -tune film <outfilename>.mp4
#
# See the H.264 encoding tutorial at
# https://trac.ffmpeg.org/wiki/Encode/H.264
# ******************************************************************************



######################### USER-DEFINED VARIABLES ###############################

paths=(
    "Frames"
)

basenames=(
    "it_"
)

extensions=(
    ".png"
)

durations=(
    0.03
)

list="ffmpeg_list.txt"

################################################################################


# Sanity check
L=${#paths[@]}
L_basenames=${#basenames[@]}
L_extensions=${#extensions[@]}
L_durations=${#durations[@]}

if [ $L_basenames != $L ] || [ $L_extensions != $L ] || [ $L_durations != $L ]
then
    echo "ERROR: mismatching dimensions of input arrays"
    exit 255
fi

if [ -f $list ]
then
    rm -rf $list
fi

for ((l=0; l<$L; ++l))
do
    path=${paths[$l]}
    path+="/"
    path+=${basenames[$l]}
    path+="*"

    for file in $path
    do
        echo "file '$file'" >> $list
        echo "duration ${durations[$l]}" >> $list
    done
done

if [ -f $list ]
then
    echo "File '$list' generated successfully"
else
    echo "Problem generating file '$list'"
    exit 255
fi
