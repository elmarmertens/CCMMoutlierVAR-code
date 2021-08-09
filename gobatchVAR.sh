#!/bin/bash
# shell script to launch matlab codes in background

wd=\'`pwd`\'
# echo $wd
preamble='cd '${wd}'; '
# echo $preamble

os=$(uname)

if [ $os = Darwin ]
then
    alias matlab='~/Applications/MATLAB_R2021a.app/bin/matlab'
    caffeinate -iw $$ &
fi


doCensorYields=true
mfile=storeVAR.m
mcmcdraws=1e3
jumpoffDate='datenum(2021,03,1)'

for modeltype in $@ 
do
    
    echo $modeltype
    # thismfile=`basename $mfile .m`

    sed '/SED-PARAMETERS-HERE/a\ 
        modeltype='\'$modeltype\'';\
        MCMCdraws='$mcmcdraws';\
	fcstNdraws= 100 * MCMCdraws;\
        jumpoffDate='$jumpoffDate';\
	doCensorYields='$doCensorYields';\
	'\
	$mfile > foobatch.m
    
    # cat foobatch.m
    matlab -nodisplay -batch "$preamble foobatch"

done

