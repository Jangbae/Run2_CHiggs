import os,sys,datetime,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from utils import *

thisDir = os.getcwd()
outputDir = thisDir+'/'

region='HT' #SR,CR --> matters only when plotting kinematics
categorize=1 #1==categorize into t/W/b/j, 0==only split into flavor
sigTrainedList=['250']#,'500','1000']

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
time='%i_%i_%i'%(cTime.hour,cTime.minute,cTime.second)

iPlotList = [#distribution name as defined in "doHists.py"

# 			'minBBdr',
# 			'aveBBdr',
# 			'deltaEta_maxBB',
# 			'FW_momentum_2',
# 			'centrality',
# 			'aveCSVpt',
# 			'HT',
# 			'minMlb',
# 			'Bjet1Pt',
# 			'mass_maxJJJpt',
# 			'MTlmet',
# 			'lepDR_minBBdr',
# 			'MET',
# #  
# 			'NPV',
			'lepPt',
			'lepEta',
			'JetEta',
			'JetPt',
# 			'NJets',
# 			'NBJets',
# 			'HTpBDT',
# 			'deltaPhi_METjets',
#			'min_deltaPhi_METjets'

# 			'HTpDNN',	
			]

isEMlist = ['E','M']
nttaglist = ['0p']
nWtaglist = ['0p']
nbtaglist = ['1','2','3p']
njetslist = ['3','4','5','6p']

if not categorize: 
	nbtaglist = ['1p']
	njetslist = ['3p']
if not categorize and 'BDT' in region: 
	nbtaglist = ['2p']
	njetslist = ['5p']

catList = list(itertools.product(isEMlist,nttaglist,nWtaglist,nbtaglist,njetslist))


count=0
for sigTrained in sigTrainedList:
	pfix='templates'
	if not categorize: pfix='kinematics_'+region
	pfix+='_M'+sigTrained+'_'+date#+'_'+time
	outDir = outputDir+pfix
	if not os.path.exists(outDir): os.system('mkdir '+outDir)
	os.chdir(outputDir)
	os.system('cp ../analyze.py doHists.py ../weights.py ../samples.py doCondorTemplates.py doCondorTemplates.sh '+outDir+'/')
	os.chdir(outDir)

	for iplot in iPlotList:
		for cat in catList:
			if skip(cat[4],cat[3]) and categorize: continue #DO YOU WANT TO HAVE THIS??
			catDir = cat[0]+'_nT'+cat[1]+'_nW'+cat[2]+'_nB'+cat[3]+'_nJ'+cat[4]
			print "Training: "+sigTrained+", iPlot: "+iplot+", cat: "+catDir
			if not os.path.exists(outDir+'/'+catDir): os.system('mkdir '+catDir)
			os.chdir(catDir)
			os.system('cp '+outputDir+'/doCondorTemplates.sh '+outDir+'/'+catDir+'/'+cat[0]+'T'+cat[1]+'W'+cat[2]+'B'+cat[3]+'J'+cat[4]+iplot+'.sh')						
	
			dict={'dir':outputDir,'iPlot':iplot,'region':region,'isCategorized':categorize,
			      'isEM':cat[0],'nttag':cat[1],'nWtag':cat[2],'nbtag':cat[3],'njets':cat[4],
			      'exeDir':outDir+'/'+catDir,'sigTrained':sigTrained}
	
			jdf=open('condor.job','w')
			jdf.write(
"""universe = vanilla
Executable = %(exeDir)s/%(isEM)sT%(nttag)sW%(nWtag)sB%(nbtag)sJ%(njets)s%(iPlot)s.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_memory = 3072
Output = condor_%(iPlot)s.out
Error = condor_%(iPlot)s.err
Log = condor_%(iPlot)s.log
Notification = Error
Arguments = %(dir)s %(iPlot)s %(region)s %(isCategorized)s %(isEM)s %(nttag)s %(nWtag)s %(nbtag)s %(njets)s %(sigTrained)s
Queue 1"""%dict)
			jdf.close()

			os.system('condor_submit condor.job')
			os.system('sleep 0.5')
			os.chdir('..')
			count+=1

print "Total jobs submitted:", count
                  
