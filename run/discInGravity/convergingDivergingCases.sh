#!/bin/bash

#Make a sequence of copies of baseCase and change one or more parameters in the
#setup files of each case

#Settings

baseCase=${1:-baseCase_6DoFRBM}

scanDir=${2:-convergingDivergingCases}

if [ ! -d "$scanDir" ];
then
    mkdir $scanDir
    cp -r $baseCase $scanDir
fi

cd $scanDir


parm1File="setCaseParms"
parm1Name=bodyDensity
parm1ShortName=rhob
parm1Range=(0.9 1.1)

parm2File="constant/dynamicMeshDict"
parm2Name=sixDoFRigidBodyMotionCoeffs.accelerationRelaxation
parm2ShortName=aRelax
parm2Range=(1 1)

parm3File="system/fvSolution"
parm3Name="PIMPLE.nOuterCorrectors"
parm3ShortName="nOC"
parm3Range=(1 1)

##########################

echo "=========================================" >> README
echo "Parameter scan generated with scanner script" >> README
echo "Date: " $(date) >> README
echo "Base case: $baseCase" >> README
echo "Parameters: $parm1Name    $parm2Name    $parm3Name" >> README
echo "=========================================" >> README
echo "Case    $parm1Name    $parm2Name    $parm3Name" >> README

for m in ${!parm1Range[*]}
do
    parm1=${parm1Range[$m]}
    parm2=${parm2Range[$m]}
    parm3=${parm3Range[$m]}
    caseName=$parm1ShortName$parm1$parm2ShortName$parm2$parm3ShortName$parm3
    if [ "${1}" = "run" ];
    then
        cd ${caseName}
        ./Allrun &
        cd ..
    elif [ ! -d "$caseName" ];
    then
        cp -r ../$baseCase $caseName
        sed -i "s/${parm1Name}=.*/${parm1Name}=${parm1}/" ${caseName}/${parm1File}
        foamDictionary -entry "${parm2Name}" -set "${parm2}" -disableFunctionEntries ${caseName}/${parm2File} 
        foamDictionary -entry "${parm3Name}" -set "${parm3}" -disableFunctionEntries ${caseName}/${parm3File} 
        echo "$caseName:    $parm1    $parm2    $parm3" >> README
    else
        echo "Directory $caseName already exists - not recreated."
    fi
done

echo "=========================================" >> README
