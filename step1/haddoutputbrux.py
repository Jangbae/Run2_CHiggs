import os,sys,datetime,time, subprocess, math
from ROOT import *

start_time = time.time()

#IO directories must be full paths
# input  = sys.argv[1]
# output = sys.argv[2]
shift = sys.argv[1]

inDir='/mnt/hadoop/store/user/jblee/CHiggs2017/'+shift+'/'
outDir='/mnt/hadoop/store/user/jblee/CHiggs2017/step1hadds/'+shift+'/'

os.system('mkdir -p '+outDir)

dirList = [
# 'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8',
# 'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8',
# 'ST_s-channel_antitop_leptonDecays_13TeV-PSweights_powheg-pythia',
# 'ST_s-channel_top_leptonDecays_13TeV-PSweights_powheg-pythia',
# 'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
# # 'SingleElectron',
# # 'SingleMuon',
# 'TTHH_TuneCP5_13TeV-madgraph-pythia8',
# 'TTTJ_TuneCP5_13TeV-madgraph-pythia8',
# 'TTTT_TuneCP5_13TeV-amcatnlo-pythia8',
# 'TTTW_TuneCP5_13TeV-madgraph-pythia8',
# 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TTToSemiLepton_HT500Njet9_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TTWH_TuneCP5_13TeV-madgraph-pythia8',
# 'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8',
# 'TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8',
# 'TTWW_TuneCP5_13TeV-madgraph-pythia8',
# 'TTWZ_TuneCP5_13TeV-madgraph-pythia8',
# 'TTZH_TuneCP5_13TeV-madgraph-pythia8',
# 'TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',
# 'TTZZ_TuneCP5_13TeV-madgraph-pythia8',
# 'TT_Mtt-1000toInf_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8',
# 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WW_TuneCP5_13TeV-pythia8',
# 'WZ_TuneCP5_13TeV-pythia8',
# 'ZZ_TuneCP5_13TeV-pythia8',
# 'ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8',
# 'ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8',

# 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToSemiLepton_HT500Njet9_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TT_Mtt-1000toInf_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8',

]

if shift=='nominalll':
	dirList.append('SingleElectron')
	dirList.append('SingleMuon')

for sample in dirList:
    outList = ['none']
    if 'TTTo' in sample: outList = ['Mtt0to700','Mtt700to1000','Mtt1000toInf']

    for outlabel in outList:

        outsample = sample+'_'+outlabel
        if outlabel == 'none': outsample = sample

        rootfiles = [x for x in os.listdir(inDir+'/'+outsample) if '.root' in x]

        print "------------ hadding Sample:",outsample,"---------------"
        print 'N root files in',outsample,'=',len(rootfiles)


        nFilesPerHadd = 900
        if 'TTToSemiLeptonic' in sample and outlabel=='Mtt0to700': nFilesPerHadd = 45
        elif 'WJetsToLNu_HT-1200To2500' in sample: nFilesPerHadd = 20
        elif 'WJetsToLNu_HT-2500ToInf' in sample: nFilesPerHadd = 13        
        onefile = inDir+'/'+outsample+'/'+rootfiles[-1]
        manyfiles = nFilesPerHadd*onefile
        lengthcheck = len('hadd -f '+outDir+'/'+outsample+'_hadd.root '+manyfiles)
        if lengthcheck > 126000.:
            toolong = lengthcheck - 126000.
            num2remove = math.ceil(toolong/len(onefile))
            nFilesPerHadd = int(nFilesPerHadd - num2remove)
            print 'Length estimate reduced from',lengthcheck,'by',toolong,'via removing',num2remove,'files for nFilesPerHadd of',nFilesPerHadd

        if len(rootfiles) < nFilesPerHadd:
            haddcommand = 'hadd -f '+outDir+'/'+outsample+'_hadd.root '
            for file in rootfiles:
                haddcommand+=' '+inDir+'/'+outsample+'/'+file
            print 'Length of hadd command =',len(haddcommand)
            subprocess.call(haddcommand,shell=True)
        else:
            for i in range(int(math.ceil(len(rootfiles)/float(nFilesPerHadd)))):
                haddcommand = 'hadd -f '+outDir+'/'+outsample+'_'+str(i+1)+'_hadd.root '

                begin=i*nFilesPerHadd
                end=begin+nFilesPerHadd
                if end > len(rootfiles): end=len(rootfiles)
                print 'begin:',begin,'end:',end

                for j in range(begin,end):
                    haddcommand+=' '+inDir+'/'+outsample+'/'+rootfiles[j]
                print 'Length of hadd command =',len(haddcommand)
                subprocess.call(haddcommand,shell=True)

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))



