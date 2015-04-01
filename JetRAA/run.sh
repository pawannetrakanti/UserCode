#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
source /osg/app/cmssoft/cms/cmsset_default.sh


cd /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/
eval `scramv1 runtime -sh`
cd -

type=$1
infile=$2
outfile=$3

#root -l -b -q /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JetRaa/CondorPbPb_Data.C+\(\"$type\",\"$infile\",\"$outfile\"\)

root -l -b -q /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JetRaa/jetmatch.C+\(\"$type\",\"$infile\",\"$outfile\"\)

# Physical path 
# /net/hidsk0001/d00/scratch/pawan/condorfiles/pbpb/Data
# /net/hidsk0001/d00/scratch/pawan/condorfiles/pbpb/ntuples

scp $outfile hidsk0001:/export/d00/scratch/pawan/condorfiles/pbpb/ntuples
rm -f $outfile

echo "Done!"

