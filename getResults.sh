#!/bin/sh
DIRECTORY="outputs/raw_outputs"
if [ ! -d "$DIRECTORY" ]; then
    mkdir $DIRECTORY
fi
cd $DIRECTORY
if [ ! -d "results" ]; then
    mkdir results
    mkdir logs
fi

cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_mid/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_102/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_106/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_110/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_114/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_118/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_122/results/*.root results/

cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_mid/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_102/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_106/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_110/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_114/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_118/results/*.root results/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_122/results/*.root results/


cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_mid/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_102/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_106/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_110/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_114/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_118/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SP_122/logs/*.logs logs/

cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_mid/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_102/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_106/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_110/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_114/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_118/logs/*.logs logs/
cp /home/w955m639/v1flow/2015_PbPb/GenerateV1/results/v1SPmc_122/logs/*.logs logs/
