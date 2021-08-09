#!/bin/bash
# shell script to launch matlab codes in background

wd=\'`pwd`\'
# echo $wd
preamble='cd '${wd}'; parpool; '
# echo $preamble

os=$(uname)

if [ $os = Darwin ]
then
    alias matlab='~/Applications/MATLAB_R2021a.app/bin/matlab'
    caffeinate -iw $$ &
fi


doCensorYields=true

for mfile in $@ 
do
    thismfile=`basename $mfile .m`

    sed '/SED-PARAMETERS-HERE/a\ 
        doCensorYields='$doCensorYields';'\
	$mfile > foobatch.m
    
    matlab -nodisplay -batch "$preamble foobatch"

done

