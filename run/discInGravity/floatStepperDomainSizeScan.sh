#!/bin/bash

#Make a sequence of copies of baseCase and change one or more parameters in the
#setup files of each case

#Settings

baseCase=${1:-baseCase_floatStepper}

scanDir=${2:-floatStepperDomainSizeScan}

if [ ! -d "$scanDir" ];
then
    mkdir $scanDir
    cp -r $baseCase $scanDir
fi

cd $scanDir


parm1File="writeMeshParmFile.py"
parm1Name=R2
parm1ShortName=R
parm1Range=(40 60 80)

parm2File="writeMeshParmFile.py"
parm2Name=nx
parm2ShortName=nx
parm2Range=(200 300 400)


##########################

echo "=========================================" >> README
echo "Parameter scan generated with scanner script" >> README
echo "Date: " $(date) >> README
echo "Base case: $baseCase" >> README
echo "Parameters: $parm1Name    $parm2Name" >> README
echo "=========================================" >> README
echo "Case    $parm1Name    $parm2Name" >> README

for m in ${!parm1Range[*]}
do
    parm1=${parm1Range[$m]}
    parm2=${parm2Range[$m]}
    caseName=$parm1ShortName$parm1$parm2ShortName$parm2
    if [ "${1}" = "run" ];
    then
        cd ${caseName}
        ./Allrun &
        cd ..
    elif [ ! -d "$caseName" ];
    then
        cp -r ../$baseCase $caseName
        sed -i "s/${parm1Name}\s*=.*/${parm1Name}=${parm1}/" ${caseName}/${parm1File}
        sed -i "s/${parm2Name}\s*=.*/${parm2Name}=${parm2}/" ${caseName}/${parm2File}
        echo "$caseName:    $parm1    $parm2" >> README
    else
        echo "Directory $caseName already exists - not recreated."
    fi
done

echo "=========================================" >> README
