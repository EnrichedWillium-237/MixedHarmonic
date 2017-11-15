#!/bin/sh
cd ..
DIRECTORY="results/v1SP_114"
if [ ! -d "$DIRECTORY" ]; then
    mkdir $DIRECTORY
    cp -r src $DIRECTORY/.
    cp GenerateV1.C $DIRECTORY/.
    cp -r EfficiencyTables $DIRECTORY/.
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_262548_262799/* > inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_262800_263230/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_263231_263359/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_263360_263379/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_263380_263614/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_tight_263615_263757/* >> inlist.dat
mv inlist.dat $DIRECTORY/.
fi
cd $DIRECTORY
if [ ! -d "results" ]; then
    mkdir results
    mkdir logs
fi
# root -l -b -q '/home/w955m639/v1flow/2015_PbPb/GenerateV1/GenerateV1.C+("v1SP_114", "inlist.dat", 50)'
root -l -b -q '/home/w955m639/v1flow/2015_PbPb/GenerateV1/GenerateV1.C+("v1SP_114", "inlist.dat", 5000)' > logs/v1SP_114.logs
cd ..
