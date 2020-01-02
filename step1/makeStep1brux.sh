#!/bin/bash

hostname
date

infilename=${1}
outfilename=${2}
inputDir=${3}
outputDir=${4}
idlist=${5}
ID=${6}
scratch=${PWD}

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_10
cd CMSSW_10_2_10
eval `scramv1 runtime -sh`
cd -

echo "setting macroDir to PWD"
macroDir=${PWD}
export PATH=$PATH:$macroDir

XRDpath=root://cmseos.fnal.gov/$inputDir

echo "Running step1 over list: ${idlist}"
for iFile in $idlist; do
    inFile=${iFile}
    if [[ $iFile == ext* ]] ;
    then
	inFile=${iFile:4}
    elif [[ $iFile == [ABCDEFWXYZ]* ]] ;
    then
	inFile=${iFile:1}
    fi

    echo "creating ${outfilename}_${iFile}.root by reading ${infilename}_${inFile}"
    root -l -b -q makeStep1.C\(\"$macroDir\",\"$XRDpath/${infilename}_${inFile}.root\",\"${outfilename}_${iFile}.root\"\)
done

echo "ROOT Files:"
ls -l *.root

# copy output to eos

NOM="nominal"
echo "cp output for condor"
#for SHIFT in nominal JECup JECdown JERup JERdown
for SHIFT in nominal
  do
  haddFile=${outfilename}_${ID}${SHIFT}_hadd.root
  hadd ${haddFile} *${SHIFT}.root
  echo "cp ${haddFile} ${outputDir//$NOM/$SHIFT}/${haddFile//${SHIFT}_hadd/}"
  cp ${haddFile} ${outputDir//$NOM/$SHIFT}/${haddFile//${SHIFT}_hadd/}
  rm *${SHIFT}.root
  rm ${haddFile}
  if [[ $haddFile == Single* ]]; then break; fi;
done

echo "done"
