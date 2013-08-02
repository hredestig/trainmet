import ghmm
import random
import GQLMixture
import GQL
import GQLCluster
import math
import sys
import numpy
import csv
import copy
import getopt
import scipy.stats
import sets

def readPathwayAnnotation(inputFile):		
  file = open(inputFile)
  pathwaysToIndividuals = {}
  individualsTopathways = {}

  for i,line in enumerate(file):
    aux = line.strip('\n')
    aux = aux.split(' ')
    #print aux
    pathway = aux[0]
    individuals = aux[1].split('|')
    #print individuals
    pathwaysToIndividuals[pathway] = individuals
    for ind in individuals:
        try:
          value = individualsTopathways[ind]
          value.append(pathway)
        except KeyError:
          individualsTopathways[ind] = [pathway]
  return pathwaysToIndividuals, individualsTopathways


def readPathwayAnnotation(inputFile):		
  file = csv.reader(open(inputFile), delimiter='\t')
  pathwaysToIndividuals = {}
  individualsTopathways = {}

  for i,line in enumerate(file):
    aux = line.strip('\n')
    aux = aux.split('\t')
    pathway = aux[0]
    individuals = aux[1].split('|')
    pathwaysToIndividuals[pathway] = individuals
    for ind in individuals:
        try:
          value = individualsTopathways[ind]
          value.append(pathway)
        except KeyError:
          individualsTopathways[ind] = [pathway]
  return pathwaysToIndividuals, individualsTopathways

 
def correlation_distribution(profileSet,profileSetTest,lag=[0]):
    logpearson = numpy.zeros((len(profileSetTest),len(profileSet)))
    lagres = numpy.zeros((len(profileSetTest),len(profileSet)))
    corr =  numpy.zeros((len(profileSetTest),len(profileSet)))
    n = len(profileSetTest[0])
    masks = []
    aux = numpy.array(range(n))
    #for l in lag:
    #    if l >= 0:
    #      masks.append((aux < n - l,aux > l - 1))
    #    else:
    #      masks.append((aux > abs(l) - 1,aux < n - abs(l)))
    print 'masks', masks
    tdist = scipy.stats.t(n-2)
    for i,s in enumerate(profileSet):
        s = numpy.array(s)
        sel = s>-9999
        s = s[sel]
        n = len(s)
        aux = numpy.array(range(n))
        masks = []
        for l in lag:
          if l >= 0:
            masks.append((aux < n - l,aux > l - 1))
          else:
            masks.append((aux > abs(l) - 1,aux < n - abs(l)))
          
        for j,s1 in enumerate(profileSetTest):
            s1 = numpy.array(s1);
            s1 = s1[sel]
            ps = []
            corrs = []
            for (mask1,mask2) in masks:
              #print len(s), len(s1), len(mask1), len(mask2), len(sel), len(s[sel & mask1]), len(s1[sel & mask2]), sum(mask1), sum(mask2),sel, mask1
              [r,p] = scipy.stats.pearsonr(s[mask1],s1[mask2]);
              ps.append(p)
              corrs.append(r)
            lagres[j,i] = lag[numpy.argmin(ps)]
            corr[j,i] = corrs[numpy.argmin(ps)]
            logpearson[j,i]=numpy.log(1-numpy.min(ps))
            
    return logpearson, lagres, corr

def zscorepearson(profileSet,profileSetTest):
    [log,lag,corr] = correlation_distrubution(profileSet,profileSetTest)
    zscore = numpy.zeros((len(corr),len(corr[0])))
    for i,cs in enumerate(corr):
      for j,c in enumerate(cs):
        zscore[i,j] = numpy.log((1+c)/(1-c))/2
    return zscore

def discretize(values):
    res = numpy.zeros((len(values),len(values[0])))
    res[values<-1.5] = 0
    res[values>=-1.5 & values < -0.5] = 1
    res[values>=-0.5 & values < 0.5] = 2
    res[values>=0.5 & values < 1.5] = 3
    res[values>=1.5] = 4
    return res

def conditionalProb(disccorr,metabolites,mapmetclasses,mettranscript):
    # Method implementing the posterior probability of correlation given the metabolite class and functional relation
    # See troyanskay for detais
    # maybe one needs pseudo counts
    
    counts = numpy.zeros((max(mapmetclasses),2,5))
    for j,met in enumreate(metabolites):
        k = mapmetclasses[met]
        for i in enumerate(len(mettranscript)):
            if mettranscript[i,j]:
              counts[k,0,disccorr[i,j]] += 1
            else:
              counts[k,1,disccorr[i,j]] += 1
                  
    summ = numpy.sum(counts,axis=2)
    for i in len(counts[0,0,:]):
        counts[:,:,i] = counts[:,:,i]/summ
    return counts
                         
def evalconditionalProb(disccorr,metabolites,mapmetclasses,conddist):
    
    res = numpy.zeros((len(disccorr),len(metabolites),2))
    for j,met in enumerate(metabolites):
        k = mapmetclasses[met]
        for i in enumerate(len(disccorr)):
          res[i,j,:] = conddist[k,:,disccorr[i,j]]
    return res
                         
def cv_samples(folds,n):
    foldSize = math.ceil(n/folds)
    rest = n%folds
    indices =range(n)
    random.shuffle(indices)
    results = []
    for i in range(folds):
      aux = []
      if i < rest:
        for j in range(foldSize+1):
           aux.append(indices.pop())
      else:
        for j in range(0,foldSize):
           aux.append(indices.pop())
      aux.sort()
      results.append(aux)
    return results

def fullModelTrain(disccorrs,metabolites,mapmetclasses,mettranscript,folds):
    # Perfom training and evaluate with cross-validations
    testfolds = cv_samples(folds,len(mettranscript))
    res = numpy.zeros((len(mettranscript),len(metabolites)))
    condfolds = []
    priorfolds = []
    #training the model
    for fold in testfold:
        train = range(len(mettranscript))
        for i in test:
            train.pop(i)
        prior = numpy.zeros((2,1))
        prior[0] = sum(sum(mettranscript[train,:]==0))
        prior[1] = sum(sum(mettranscript[train,:]==1))
        prior = prior/sum(prior)
        priorfolds.append(prior)
        cprobfold = []
        for disccorr in disccorrs:
            cprobfold.append(conditionalProb(disccorr[train,:],metabolites,mapmetclasses,mettranscript[train,:]))
        condfolds.append(folds)
    #evaluating the model
    res = numpy.ones((len(disccorr),len(metabolites),2))
    for i,fold in enumerate(testfold):
        for t in fold:
            for j,m in enumerate(metabolites):
              k = mapmetclasses[met]
              for i in enumerate(len(disccorrs)):
                res[t,j,:]*=conddist[k,:,disccorr[t,j]]
              res[t,j,0]*=priorfolds[i][0]
              res[t,j,1]*=priorfolds[i][1]
              res[t,j,:]=res[t,j,:]/sum(res[t,j,:])
    return res

def bradley_method(met_files,genes_files,path_met_files,path_genes_files,class_file,folds):
  
  [classes2metabolites, metabolites2classes]=readPathwayAnnotation(class_file)
  [pathways2metabolites, metabolites2pathways]=readPathwayAnnotation(path_met_files)
  [pathways2genes, genes2pathways]=readPathwayAnnotation(path_met_files)

  transcriptMetabolite = numpy.zeros((len(genes2pathways),len(metabolites2classes)))

  for i,g in enumerate(genes2pathways.ids()):
    geneSet = sets.Set(genes2pathways[g])
    for j,m in enumerate(metabolites2pathways.ids()):
      if len(geneSet.intersection(metabolites2pathways[m])) > 0:
        transcriptMetabolite[i,j]=1
      else:
        transcriptMetabolite[i,j]=0   

  disccorrs = []
  metabolites = []
  
  for i,f in enumerate(met_files):
    
    profileSetTest = GQL.ProfileSet()
    profileSetTest.ReadDataFromCaged(genes_files[i],end=1)

    profileSet = GQL.ProfileSet()
    profileSet.ReadDataFromCaged(f,end=1)

    metabolites.append(profileSet.acc)
    
    values = zscorepearson(profileSet,profileSetTest)
    values = discretize(values)   
    disccorrs.append(values)


  res = fullModelTrain(disccorrs,metabolites,metabolites2classes,transcriptMetabolite,folds)  
                                   

def estimate_models(profileSet,outFile,cyclic,states):
    profiles = []
    models = []

    for l in range(len(profileSet.ghmm_seqs)):
        profileSet.ghmm_seqs.setWeight(l,0.0)

    viterbi = []
        
    for l in range(len(profileSet.ghmm_seqs)):
          [currentProfile,model,vit]=train_hmm(profileSet,states,cyclic, l)
          profiles.append(currentProfile)
          models.append(model)
          viterbi.append(vit)
          #print 'original', vit
    return [profiles,models,profileSet,viterbi]

def estimate_models_mixture(profileSet,outFile,cyclic):
    
    profiles = [] 
    models = []

    for l in range(len(profileSet.ghmm_seqs)):
        profileSet.ghmm_seqs.setWeight(l,0.0)
        
    viterbi = []
        
    for l in range(len(profileSet.ghmm_seqs)):
        modelsMix = []
        viterbiMix = []
        profilesMix = []
        states = range(2,len(profileSet.ghmm_seqs[0])/2+1)
        if cyclic:
          states.append(2)
          states.append(3)
        cyclicmodels = numpy.zeros((len(states),1));
        if cyclic:
          cyclicmodels[-1] = 1
          cyclicmodels[-2] = 1
        for i in range(len(states)):
           [currentProfile,model,vit]=train_hmm(profileSet,states[i],cyclicmodels[0], l)
           profilesMix.append(currentProfile)
           modelsMix.append(model)
           viterbiMix.append(vit)
        models.append(modelsMix)
        viterbi.append(viterbiMix)
        profiles.append(profilesMix)
    return [profiles,models,profileSet,viterbi]



def train_hmm(profileSet,states,cyclic, l):

    currentProfile = GQLCluster.ProfileClustering()
    currentProfile.setProfileSet(profileSet)
          #profiles.append(currentProfile)
          
    param = random_models_args(states)
    currentProfile.randomModels(param,cyclic=cyclic,noise=0.05)

    #set weight for training only with observation i
    profileSet.ghmm_seqs.setWeight(l,1.0)
    #soft train model
    model = currentProfile.models[0]
    model.baumWelch(profileSet.ghmm_seqs, 5, 0.1)      
    profileSet.ghmm_seqs.setWeight(l,0.0)
    model.write('models/'+outFile+'_'+profileSet.acc[l]+'.xml')
    #models.append(model)

    vit = model.viterbi(profileSet.ghmm_seqs)[0]
    #viterbi.append(vit[l])
    #print 'antes vit', vit
    return [currentProfile,model,vit[l]] 

def estimate_models_pathways(fileList,outFile,cyclic):
    pathways = []
    profiles = []
    models = []
    fileList = open(fileList,'r')
    for l in fileList:
        file = l.strip('\n')
        fileAux = file.split('/')
        pathway_name = fileAux[-1]
        pathways.append(fileAux[-1])
        profileSet = GQL.ProfileSet()
        profileSet.ReadDataFromCaged(file,end=1)
        
        currentProfile = GQLCluster.ProfileClustering()
        currentProfile.setProfileSet(profileSet)
        profiles.append(currentProfile)
          
        param = random_models_args()
        currentProfile.randomModels(param,cyclic=cyclic,noise=0.05)
        #soft train model
        model = currentProfile.models[0]
        model.baumWelch(profileSet.ghmm_seqs, 5, 0.1)      
        model.write('models/'+outFile+'_'+pathway_name+'.xml')
        models.append(model)
          
    return [profiles,models,pathways]

def diagonal(diag):
    values = numpy.zeros((len(diag),len(diag)))
    for i in range(len(values)):
        values[i,i] = diag[i]
    return values

def estimateLikelihoods(profileSet,profileSetTest,profileSetTestNegative,models,viterbiref=[],mixture=0,pearsonmix=0):
      
    lls = []
    lags = []
    vits = []
    lagns = []
    llns = []

    if mixture:
      llaux = numpy.zeros((len(profileSetTest),2*len(models[0])))
    if pearsonmix:
      llaux = numpy.zeros((len(profileSetTest),2*len(models[0])+1))  
      [corr,lagcorr,corr] = numpy.exp(correlation_distribution(profileSet,profileSetTest))
    

    for l,model in enumerate(models):
      #print likelihood of all sequences (should be maximal for sequence i with model i)

      if mixture:
        for i,m in enumerate(model):
          llaux[:,i] = numpy.exp(m.loglikelihoods(profileSetTest.ghmm_seqs))
        for m in model:
          i = i + 1;
          llaux[:,i] = numpy.exp(m.loglikelihoods(profileSetTestNegative.ghmm_seqs))

          sumcols = numpy.sum(llaux,axis=0)
          if pearsonmix:
            llaux[:,-1] = corr[:,l]
            sumcols[-1] = numpy.sum(sumcols[:-2]) # pearson has half of the weight
          sumcols = sumcols/sum(sumcols)            

        for s in sumcols:
            llaux[:,i] = llaux[:,i]*s;

        ll = numpy.log(numpy.sum(llaux,axis=1))
        indices = numpy.argmax(llaux,axis=1)
        lln = ll

        print 'viterbi'
        viterbi = [0]*len(profileSetTestNegative);
        for j,i in enumerate(indices):
          #print 'best index', i
          if i == 2*len(model):
            viterbi[j] = [1]*len(profileSetTestNegative.ghmm_seqs[j])
          elif i >= len(model): # negative model
            #print i, len(model)
            modelaux = model[i-len(model)]  
            viterbi[j] = modelaux.viterbi(profileSetTestNegative.ghmm_seqs[j])[0]
          else: # positive model
            modelaux = model[i]  
            viterbi[j] = modelaux.viterbi(profileSetTest.ghmm_seqs[j])[0]
        viterbin = viterbi
                            
                            
      else:
        ll = model.loglikelihoods(profileSetTest.ghmm_seqs)
        lln = model.loglikelihoods(profileSetTestNegative.ghmm_seqs)
        viterbi = model.viterbi(profileSetTest.ghmm_seqs)[0]
        viterbin = model.viterbi(profileSetTestNegative.ghmm_seqs)[0]
      finalviterbi = [];
      lag = [];
      lagn = [];
      if len(viterbi) > 0:
        viterbireference = viterbiref[l]
        try:
          len(viterbireference[0])
        except TypeError:
          viterbireference = [viterbireference]
        lagreference = []
        for i in viterbireference:
          lagreference.append(change(i,0,1));
        #print 'ref', viterbireference, lagreference
      else:
        lagreference = 0

      for i in range(len(viterbi)):
        try:
          finalviterbi.append(viterbi[i])
          if mixture:
            index = indices[i]
          else:
            index = 0
          lag.append(change(viterbi[i],0,1)-lagreference[index])
          #print 'final', viterbi[i], change(viterbi[i],0,1)-lagreference[index]
          #print 'mean diff', numpy.mean(numpy.array(viterbi[i])-numpy.array(viterbireference))
          finalviterbi.append(viterbin[i])
          lagn.append(change(viterbin[i],0,1)-lagreference[index])
        except TypeError:
          print viterbi[i], viterbin[i]
          lag.append(0)
          lagn.append(0)
          finalviterbi.append([])
                        
      lls.append(ll)
      llns.append(lln)
      lags.append(lag)
      lagns.append(lagn)
      vits.append(finalviterbi)
      print "best sequence for model ",str(l), "is",str(numpy.argmax(ll))

    return [lls,llns,lags,lagns,vits]
      

def change(viterbi,first,second):
    for i,v in enumerate(viterbi):
        if i==(len(viterbi)-1):
          return 0;
        if(v == first) & (viterbi[i+1] == second):
          return i

def argmax(iterable):
    return max((x,i) for i,x in enumerate(iterable))[1]

class random_models_args:
	def __init__(self, states):
		self.number_of_states = '['+str(states)+']'
		self.number_of_models = 1
		self.total_duration = 8
		self.parameters = 0
		self.default_variance = 3.0
		self.default_mean = 4.0
		self.noiseModel = 0
                self.dimension = dim

def readTfAnnotation(inputFile):		
  file = csv.reader(open(inputFile), delimiter='\t')

  annot = []
  for i,line in enumerate(file):
    if i == 0:
      tfs = line	    
    else:
      annot.append([int(v) for v in line[1:]])   
  annot = numpy.array(annot)
  print len(annot)
  return tfs, annot # numpy.transpose(annot)

def results(outFile,ids,profileSetTest,lls,llns,lgs,lgns,vits):
    out = open(outFile+'.res','w');
    header = 'Gene'
    for i in ids:
      header = header + '\t' + i
    out.write(header+'\n')
    values = numpy.array(lls)
    values = numpy.transpose(values)

    for j,line in enumerate(values):
      out.write(profileSetTest.genename[j]+'\t'+('\t'.join(["%+.4f" %i
 for i in line]))+'\n')

    out = open(outFile+'_n.res','w');
    header = 'Gene'
    for i in ids:
      header = header + '\t' + i
    out.write(header+'\n')
    values = numpy.array(llns)
    values = numpy.transpose(values)

    for j,line in enumerate(values):
      out.write(profileSetTest.genename[j]+'\t'+('\t'.join(["%+.4f" %i
 for i in line]))+'\n')
      
    out = open(outFile+'.lag','w');
    header = 'Gene'
    for i in ids:
      header = header + '\t' + i
    out.write(header+'\n')
    values = numpy.array(lags)
    values = numpy.transpose(values)
    #print len(values), len(values[0]), len(lls), len(lls[0])
    for j,line in enumerate(values):
      out.write(profileSetTest.genename[j]+'\t'+('\t'.join([str(i)
 for i in line]))+'\n')
      
    out = open(outFile+'_n.lag','w');
    header = 'Gene'
    for i in ids:
      header = header + '\t' + i
    out.write(header+'\n')
    values = numpy.array(lagns)
    values = numpy.transpose(values)
    for j,line in enumerate(values):
      out.write(profileSetTest.genename[j]+'\t'+('\t'.join([str(i)
 for i in line]))+'\n')
            
    #out = open(outFile+'.vit','w');
    #header = 'Gene'
    #for i in ids:
    #  header = header + '\t' + i
    #out.write(header+'\n')
    #values = numpy.array(vits)
    #values = numpy.transpose(values,axes=[1,0,2])
    #names = ['met']
    #names = profileSetTest.genename
    #for j,line in enumerate(values):
    #  out.write(names[j]+'\t'+('\t'.join([str(i) for i in line]))+'\n')   
    
  

import numpy

def usage():
	"""
	USage information
	"""
	print """
Usage: train_rank.py [options] <input_data_file_metabolites.txt> <input_data_file_expression.txt> <prefix_out_file>
Examples:

Options:
	-h, --help			Print this help message
	-c, --cyclic		        Use cyclic model
        -s, --states                    Number of states on HMM (default =3)
        -m, --model                     Predefine models <model_files>
        -p, --pathways                  Build models from pathway specific metaabolites
        -x, --mixture                   Build mixture of HMMs
        --pearson                       Probabilistic Pearson
        --mixpearson                    Mixture with a HMMs and pearson


"""


if __name__ == '__main__':

    states = 3
    dim = 1
    cyclic = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h:c:s:m:p:x', ['help','bradley','cyclic','states','model','pathways','mixture','pearson','pearsonlag','mixpearson'])
    except getopt.error, err:		
        usage()
        sys.exit(2)
                    
    fileNameQueries = sys.argv[-3]
    fileName = sys.argv[-2]
    outFile = sys.argv[-1]
    defautLags = [-4,-3,-2,-1,0,1,2,3,4]
    mixture = 0
    pearsonmix = 0

    type = 'N'
                                  
    for o, a in opts:
            if o in ('-h', '--help'):
                    usage()
                    sys.exit(0)
            elif o in ('-c', '--cyclic'):
                    cyclic = 1
            elif o in ('-s', '--states'):
                    states = int(a)
            elif o in ('-x', '--mixture'):
                    type = 'x'
                    mixture = 1
                    states = 0
            elif o in ('-p', '--pathway'):
                    type = 'p'
            elif o in ('--pearson'):
                    type = 'pc'
                    pearsonlags = [0]
            elif o in ('--pearsonlag'):
                    type = 'pc'
                    pearsonlags = defautLags
            elif o in ('--mixpearson'):
                    type = 'x'
                    pearsonmix = 1
                    mixture = 1
            elif o in ('-m', '--model'):
                    model = a
                    type = 'm'
            elif o in ('--bradley'):
                    type = 'b'
                    numberConditions = int(a)

    if type in 'b':
      # read bradley method
      met_files =  sys.argv[2:2+numberConditions]
      gene_files = sys.argv[2+numberConditions:2+2*numberConditions]
      path_met_files = sys.argv[-3]
      path_genes_files = sys.argv[-2]
      class_file = sys.argv[-1]      
      res = bradley_method(met_files,genes_files,path_met_files,path_genes_files,class_file,folds)

    else:
                    
      print 'type', type  
      #read data
      profileSetTest = GQL.ProfileSet()
      profileSetTest.ReadDataFromCaged(fileName,end=1)

      #read data
      profileSetTestNegative = GQL.ProfileSet()
      profileSetTestNegative.ReadDataFromCaged(fileName,end=1,negative=1)

      profileSet = GQL.ProfileSet()
      profileSet.ReadDataFromCaged(fileNameQueries,end=1)

      viterbi = []

      # set sequence weights according to posterior
      if type in 'p':
        [profiles,models,ids] = estimate_models_pathways(profileSet,outFile,cyclic)      
      elif type in 'm':
        profileClustering = GQLCluster.ProfileClustering()
        profileClustering.setProfileSet(profileSet)
        profileClustering.readModels(models,'xml')
        models = profileClustering.models
        ids = []
        for i in len(models):
          ids.append('Model '+str(i))
      elif type in 'x':
         [profiles,models,profileSet,viterbi] = estimate_models_mixture(profileSet,outFile,cyclic)
      elif 'pc' in type:
          pass
      else:
         [profiles,models,profileSet,viterbi] = estimate_models(profileSet,outFile,cyclic,states)
      ids = profileSet.genename
      if 'pc' in type:
        [lls,lags,corr] =  correlation_distribution(profileSet,profileSetTest,pearsonlags)
        lls = numpy.transpose(lls)      
        llns = lls
        lags =  numpy.transpose(lags)      
        lagns = lags
        vits = [] 
      else:
        [lls,llns,lags,lagns,vits] = estimateLikelihoods(profileSet,profileSetTest,profileSetTestNegative,models,viterbiref=viterbi,mixture=mixture,pearsonmix=pearsonmix)

      results(outFile,ids,profileSetTest,lls,llns,lags,lagns,vits)
      
