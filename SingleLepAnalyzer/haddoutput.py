import os,sys,datetime,time
from ROOT import *

start_time = time.time()

#IO directories must be full paths
shift = sys.argv[1]

inputDir = '/user_data/jlee/chargedHiggs/Analysis_2017_2018/Step1Samples/FWLJMET102X_1lep2017_052219_step1/'
outputDir = '/user_data/jlee/chargedHiggs/Analysis_2017_2018/Step1Samples/FWLJMET102X_1lep2017_052219_step1_hadd/'+shift+'/'

os.system('mkdir -p '+outputDir)

dirList = [
# 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TT_Mtt-1000toInf_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_s-channel_antitop_leptonDecays_13TeV-PSweights_powheg-pythia',
# 'ST_s-channel_top_leptonDecays_13TeV-PSweights_powheg-pythia',
# 'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8',
# 'TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8',
# 'TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
# 'TTTT_TuneCP5_13TeV-amcatnlo-pythia8',
# 'ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8',
# 'ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8',
# 
# 'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
# 
# 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 
# 'WW_TuneCP5_13TeV-pythia8',
# 'WZ_TuneCP5_13TeV-pythia8',
# 'ZZ_TuneCP5_13TeV-pythia8',
# 
# 'QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8'
    ]

if shift == 'nominal':
#     These don't need to be run for shifted directories
    dirList.append('SingleElectron')
    dirList.append('SingleMuon')

          
for sample in dirList:
    if sample == 'SingleElectron': continue
    print inputDir+'/'+sample
    rootfiles = [x for x in os.listdir(inputDir+'/'+sample+'/'+shift+'/') if '.root' in x]
    print '##########'*15
    print 'HADDING:', sample,len(rootfiles)
    print '##########'*15
    nFilesPerHadd = 999
    if sample=='TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8': nFilesPerHadd = 20
    if sample=='TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8': nFilesPerHadd = 20

    if len(rootfiles) < nFilesPerHadd:
        haddcommand = 'hadd -f '+outputDir+'/'+sample+'_hadd.root '
        for file in rootfiles:
            haddcommand+=' '+inputDir+'/'+sample+'/'+shift+'/'+file
        print haddcommand
        os.system(haddcommand)
    else:
        for i in range(int(len(rootfiles)/nFilesPerHadd)+1):
            haddcommand = 'hadd -f '+outputDir+'/'+sample+'_'+str(i+1)+'_hadd.root '
            begin=i*nFilesPerHadd
            end=begin+nFilesPerHadd
            if len(rootfiles) < end: end=len(rootfiles)
            for j in range(begin,end):
                haddcommand+=' '+inputDir+'/'+sample+'/'+shift+'/'+rootfiles[j]
            print haddcommand
            os.system(haddcommand)

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))



