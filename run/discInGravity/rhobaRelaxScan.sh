#!/bin/bash

#Make a sequence of copies of baseCase and change one or more parameters in the
#setup files of each case

#Settings

baseCase=${1:-baseCase_6DoFRBM}

scanDir=${2:-rhobaRelaxScan}

if [ ! -d "$scanDir" ];
then
    mkdir $scanDir
    cp -r $baseCase $scanDir
fi

cd $scanDir

#parm1File="initUcase/writeMeshParmFile.py"
#parm2File="initUcase/writeMeshParmFile.py"
#parm1Name="nx"
#parm2Name="ny"
#parm1Range=(200 200 200 200 200)
#parm2Range=(60 60 60 60 60)

parm3File="setCaseParms"
parm3Name=bodyDensity
parm3ShortName=rhob
parm3Range=(0.1 0.5 0.9)
#parm3Range=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1)

parm4File="constant/dynamicMeshDict"
parm4Name=accelerationRelaxation
parm4ShortName=aRelax
parm4Range=(0.5 0.9)
#parm4Range=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

##########################

echo "=========================================" >> README
echo "Parameter scan generated with aRelaxNOuterScanner" >> README
echo "Date: " $(date) >> README
echo "Base case: $baseCase" >> README
echo "Parameters: $parm3Name    $parm4Name" >> README
echo "=========================================" >> README
echo "Case    $parm3Name    $parm4Name" >> README

for m in ${!parm3Range[*]}
do
    parm3=${parm3Range[$m]}

    for n in ${!parm4Range[*]}
    do
#        parm1=${parm1Range[$n]}
#        parm2=${parm2Range[$n]}
        parm4=${parm4Range[$n]}
#        caseName=$(printf "%03d\n" "$n")
        caseName=$parm3ShortName$parm3$parm4ShortName$parm4
        if [ "${1}" = "run" ];
        then
            cd ${caseName}
            ./Allrun &
            cd ..
        elif [ ! -d "$caseName" ];
        then
            cp -r ../$baseCase $caseName
#            sed -i "s/${parm1Name}\s*=.*/${parm1Name}=${parm1}/" ${caseName}/${parm1File}
#            sed -i "s/${parm2Name}\s*=.*/${parm2Name}=${parm2}/" ${caseName}/${parm2File}
            sed -i "s/${parm3Name}\s*=.*/${parm3Name}=${parm3}/" ${caseName}/${parm3File}
#            sed -i "s/${parm3Name}\s.*;/${parm3Name}\t${parm3};/" ${caseName}/${parm3File}
            sed -i "s/${parm4Name}\s.*;/${parm4Name}\t${parm4};/" ${caseName}/${parm4File}
            echo "$caseName:    $parm3    $parm4" >> README
        else
            echo "Directory $caseName already exists - not recreated."
        fi
    done
done

echo "=========================================" >> README
