from labblouin import PDBnet,FASTAnet
import numpy as np
import pymol
import itertools
from igraph import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools
from scipy.stats import gaussian_kde
from collections import OrderedDict
import csv
import random
import scipy.stats as stats
from random import gauss
import time
import operator

#Author: Qazi Zaahirah Dalhousie University

class PDBChains:
    def __init__(self):
        filenames = None
    def ReadPDBfromTxt_GetallPDBChains(self,txtInput,inputFolder,OutputFolder):
        """This function gets all the chains of the pdb files which are read from a text file and copies it 
        in a folder"""
        mylist = []
        with open(txtInput) as textFile:
            lines = textFile.readlines()
            for line in lines:
                mylist.append(line[:4])
                
        for x in mylist: 
            p = PDBnet.PDBstructure(inputFolder+x+'.pdb')
            chainNames = p.GetChainNames()
            for ch in chainNames:
                newName = x+ch+'.pdb'
                print newName
                c = p.GetChain(ch)
                c.WriteAsPDB(OutputFolder+newName)                
            
    def GetallPDBChains(self,inputFolder,outputFolder):
        os.chdir(inputFolder)
        fileNames = os.listdir('.')
        count =1
        for f in fileNames:
            print f
            p = PDBnet.PDBstructure(f)
            chainNames= p.GetChainNames()
            print chainNames
            loop = len(chainNames)
            for i in range(0,loop):
                newName = f[:-4]+chainNames[i]+'.pdb'
                print newName
                c = p.GetChain(chainNames[i])
                c.WriteAsPDB(outputFolder+newName)        
        
class ResCent:
    """
    This class has the functions that calculates the centrality of the FCM
    """
    def __init__(self, pdbfile=None):

        if pdbfile == None:
            self.pdb = None
        else:
            self.pdb = pdbfile

        self.cent   = None 
    
    def calculate_centrality(self, matrix,centType,aggregate):
    
        '''
        Calculates the centrality values for a PDB structural set from its frequency contact matrix.
    
        :param matrix: a numpy matrix -the frequency contact matrix, calculated using PDBnet.PDBstructure.FreqCM()
        :return: a list of centrality values, with list length = protein length (in residues).
                 0 values are initialized in this list originally, then changed
                 if that residue has a centrality value above 0.
                 list order = residue order.
        '''
        shape = matrix.shape[0]
        print shape
    
        combos = [combo for combo in itertools.combinations(range(shape), 2)]
        # make a list of tuples, with each element being:
        # (resA, resB, freq)
        # if the freq is between 0 and 1
        
        if aggregate==True:
            include_combos = [(c[0],c[1],matrix[c[0]][c[1]]) for c in combos if 0 <= matrix[c[0]][c[1]] <= 1]
            #include_combos = [(c[0],c[1],(1-math.fabs(((matrix[c[0]][c[1]])-0.5)/0.5))) for c in combos if 0 < matrix[c[0]][c[1]] <= 1]
            #print include_combos
        else:    
            include_combos = [(c[0],c[1],matrix[c[0]][c[1]]) for c in combos if matrix[c[0]][c[1]]==1]
            
        g = Graph()
        ids = {}
    
        id = 0
    
        for data in include_combos:
    
            # Add vertex 1 (if it doesnt already exist)
            if data[0] not in ids:
                g.add_vertex(name=data[0])
                ids[data[0]] = id
                id += 1
            # Add vertex 2 (if it doesnt already exist)
            if data[1] not in ids:
                g.add_vertex(name=data[1])
                ids[data[1]] = id
                id += 1
    
            # Add the edge between these two vertices
            g.add_edge(ids[data[0]],ids[data[1]],weight=data[2])
    
        # return the eigenvector centrality value (normalized (0,1)) for each vertex (residue)
        if centType=='eigen':
            evcent = g.evcent(scale=True, weights=g.es['weight'])
        elif centType=='degree':
            evcent = g.degree(type="out")
    
        vx = VertexSeq(g)
    
        # vl list is the actual vertex name (which corresponds to the residue in sequence)
        vl = [vertex.attributes()['name'] for vertex in vx]
    
        vertexinfo = zip(vl, evcent) # (residue index, centrality value)
    
        centralitylist = [0.00] * matrix.shape[0] # list of zeros, length of the protein
    
        # iterate through the vertexinfo list. add to the centralitylist (at the proper index)
        # the centrality value for that residue. Remember, some values will remain 0, because
        # they were not nodes in the graph.
        for i in vertexinfo:
            centralitylist[i[0]] = i[1]
    
        #for i in range(len(centralitylist)):
            #g.vs[i]['color'] = [centralitylist[i], 0, 0]
    
        g.vs['label'] = g.vs['name']
    
        self.cent = centralitylist
        return centralitylist
    def calc_cent_onlyforhomo(self, matrix,centType):
        
            '''
            This function is a cleaner version of calculate_centrality
            Param:
            matrix: is the frequency contact matrix
            centType: is the type of centrality (there are two options 'eigen' 'degree' )
            '''
            shape = matrix.shape[0]
        
            g = Graph()
    
            for i in range(shape):
                g.add_vertex(res=i)
    
            combo_tuples = [combo for combo in itertools.combinations(range(shape), 2)]
    
            for combo in combo_tuples:
    
                #value = matrix[combo[0]][combo[1]]
                value=1-(math.fabs(((matrix[combo[0]][combo[1]])-0.5)/0.5))
    
                if value > 0 and value < 1:
    
                    g.add_edge(combo[0], combo[1], weight=value)
    
            # return the eigenvector centrality value (normalized (0,1)) for each vertex (residue)
            if centType=='eigen':
                evcent = g.evcent(scale=True, weights=g.es['weight'])
            elif centType=='degree':
                evcent = g.degree(type="out")
    
            # round the centrality value to two decimal places.
            centralitylist = [float('%.2f' % x) for x in evcent]
    
            self.cent = centralitylist
            return centralitylist    
    def graphcomponent(self,matrix):
        """
        This fucntion is used to get the components in the graph.
        Param: 
        Matrix: frequency contact matrix 
        """
        print 'this'
        d={}
        shape = matrix.shape[0]
        combos = [combo for combo in itertools.combinations(range(shape), 2)]
        includecombos=[]
        for i,values in enumerate(matrix):
            mylist=[]            
            for j,k in enumerate(values):
                if k>0:
                    mylist.append(j)
            d.update({i:mylist})    
        return d           
    def getRoots(self,aNeigh):
        def findRoot(aNode,aRoot):
            while aNode != aRoot[aNode][0]:
                aNode = aRoot[aNode][0]
            return (aNode,aRoot[aNode][1])
        myRoot = {} 
        for myNode in aNeigh.keys():
            myRoot[myNode] = (myNode,0)  
        for myI in aNeigh: 
            for myJ in aNeigh[myI]: 
                (myRoot_myI,myDepthMyI) = findRoot(myI,myRoot) 
                (myRoot_myJ,myDepthMyJ) = findRoot(myJ,myRoot) 
                if myRoot_myI != myRoot_myJ: 
                    myMin = myRoot_myI
                    myMax = myRoot_myJ 
                    if  myDepthMyI > myDepthMyJ: 
                        myMin = myRoot_myJ
                        myMax = myRoot_myI
                    myRoot[myMax] = (myMax,max(myRoot[myMin][1]+1,myRoot[myMax][1]))
                    myRoot[myMin] = (myRoot[myMax][0],-1) 
        myToRet = {}
        for myI in aNeigh: 
            if myRoot[myI][0] == myI:
                myToRet[myI] = []
        for myI in aNeigh: 
            myToRet[findRoot(myI,myRoot)[0]].append(myI) 
        return myToRet     
       
        
    def drawHistogram(self, m):
        """This is a simple function that takes a matrix as an input and
        flatens the matrix before plotting a histogram"""
        #b= np.arange(0.1,1,0.001) #make a list from 0.1 -1 with an increment of 0.001
        #change into 1D array 
        m2= m.ravel()
        #copy all the elements except 0.0 and 1.0
        m3 = m2[m2 != 0.0]
        m4= m3[m3!=1.0]
        m3=filter(lambda a: a != 2, m2)
        plt.hist(m4)
        plt.title('weight histogram')
        plt.xlabel('weight values')
        plt.ylabel('frequency')
        plt.show()
    def normalizeVector(v):
        mx = max(v)
        mn = min(v)
        centrality = []        
        for j in v:
            
            x = (j-mn)/(mx-mn)
            centrality.append(x)
            
        return centrality       
        
    

class ResEntropy:
    
    def __init__(self, pdbfile=None):

        if pdbfile == None:
            self.pdb = None
        else:
            self.pdb = pdbfile
    
        self.cent   = None 
    def calcContactEntropy (self, matrix):
        """This function just takes a matrix as an input and calculates the entropy of each row which corresponds to a residue """
        
        entropylist = []
        for i in range(matrix.shape[0]):
            
            entropylist.append(-(sum([matrix[i][j]*math.log(matrix[i][j],2) for j in range(matrix.shape[1]) if matrix[i][j] > 0])))
            
        return entropylist        
    def GetEntropyForAA(self,aalis,homolis,entropy,chain,compPath):
        """This Function is used to calculate the entropy values of 20 amino acids
        Param:
        aalis: the list of amino acids and has a length of 20
        homolis: the index of the residues that are homologous (structurally) in the entire dataset
        entropy:  a vector that contains the entropy values of the homologous residues. The length of homolis and entropy is the same
        chain: it is the chain of the PDB whose entropy for amino acids needs to be calculated
        compPath: is the path for writing the dictionary in a text file """
        newchain=[]
        for i,j in enumerate(chain):
            
            for k in homolis:
                
                if i==k:
                    newchain.append(j)
        newlist=[]
        for i,j in enumerate(entropy):
            
            newlist.append((newchain[i],j))
             
        #newdict={aa:[] for aa in aalis}
        newdict=OrderedDict()
        for aa in aalis:
            newdict.update({aa:[]})        
        
        for aa, entropy in newlist:
            
            newdict[aa].append(entropy)
            
        filepath =compPath+'.txt'
        f = open(filepath,'w')
        f.write(str(newdict))
        return newdict
        
    def BoxPlotEntropyForAA(self,aalis,homolis,entropy,chain,compPath,title):
        """This function takes entropy as a input creates a dictionary in which key correspond to 20 amino acids and the values are
        the respective values of the amino acids and then generates a box plot for the values of each amino acid"""
        newchain=[]
        for i,j in enumerate(chain):
            
            for k in homolis:
                
                if i==k:
                    newchain.append(j)
        newlist=[]
        
        for i,j in enumerate(entropy):
            
            newlist.append((newchain[i],j))
        newdict=OrderedDict()
        for aa in aalis:
            newdict.update({aa:[]})
        print newdict
        for aa, entropy in newlist:
            
            newdict[aa].append(entropy)
        labels=newdict.keys()
        print newdict
        data=[]
        
        for res in newdict.keys():
            
            mean=np.mean(newdict[res])
            sd=np.std(newdict[res])
            data1=newdict[res]
            data.append(data1)
        yticks = np.arange(0,11,2)        
        fig = plt.figure(1, figsize=(30, 20))
        fig.suptitle(title,fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel("Amino Acids", fontsize= 20)
        ax.set_ylabel("Entropy Values", fontsize =20)
        ax.set_yticks(yticks)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        bp = ax.boxplot(data,labels=labels)
        fig.savefig(compPath+'.png', bbox_inches='tight')
        filepath =compPath+'.txt'
        f = open(filepath,'w')
        f.write(str(newdict)) 
        return newdict
    def BoxPlotAggregateContactEntropyfromtxt(self,path,dictlist):
        """ This function takes dictionaries as input aggregates the dictionaries and forms the box plot"""
        f = open(path+dictlist[0]+'.txt','r')
        my_dict0 = eval(f.read())
        f = open(path+dictlist[1]+'.txt','r')
        my_dict1 = eval(f.read())
        f = open(path+dictlist[2]+'.txt','r')
        my_dict2 = eval(f.read())
        f = open(path+dictlist[3]+'.txt','r')
        my_dict3 = eval(f.read())        
        
        aggdict={a:[] for a in my_dict3.keys()}
        print aggdict
        for keys, value in my_dict0.iteritems():
            aggdict[keys].extend(value)
        for keys, value in my_dict1.iteritems():
            aggdict[keys].extend(value)           
        for keys, value in my_dict2.iteritems():
            aggdict[keys].extend(value)       
        for keys, value in my_dict3.iteritems():
            aggdict[keys].extend(value)                    
        
        data=[]        
        labels=aggdict.keys()           
        for res in aggdict.keys():
                   
                mean=np.mean(aggdict[res])
                sd=np.std(aggdict[res])
                #data1=[mean,sd]
                data1=aggdict[res]
                data.append(data1)
                

        plt.boxplot(data)
        plt.boxplot(data,labels=labels)
        plt.ylim(-1,10)
        plt.show()
        fig = plt.figure(1, figsize=(20, 15))
        ax = fig.add_subplot(111)
        bp = ax.boxplot(data,labels=labels)
        fig.savefig(path+'aggregate.png')
    def GetMeanEntropy_AA(self,AAdict,aalis):
        """Get the mean entropy for all amino acids and return the dictionary with keys as amino acids
           and the value is the mean"""
        meanAAlis = {x:0 for x in aalis}
        for key, values in AAdict.iteritems():
            mymean = mean(values)
            meanAAlis[key]= mymean
        return meanAAlis
    def rank_Plot_MeanEntropyForAA(self,newdict,aalis,plotPath,filename,count):
        """Input : Dictionary
        Output: ranking of the amino acids on the basis of the mean of entropy. This also Plots the mean entropy of amino acids"""
        meanAAlis = {x:0 for x in aalis}
        for key, values in newdict.iteritems():
            mymean = mean(values)
            meanAAlis[key]= mymean
            
        sorted_meanAAlis = sorted(meanAAlis.items(),key = operator.itemgetter(1))
        data = [x[1] for x in sorted_meanAAlis]
        mylabels = [x[0] for x in sorted_meanAAlis]
        
        xaxis = np.arange(0,20,1)
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle('Ranking of Mean Entropy for amino acids'+filename,fontsize= 30)
        ax = fig.add_subplot(111)
        #this is important so that the ticks are places at correct distance
        ax.set_xticks(xaxis)
        ax.set_ylabel('Mean of Entropy values', fontsize= 20)
        ax.set_xlabel('Amino Acids', fontsize =20)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        ax.set_xticklabels(mylabels)
        ax.bar(xaxis,data,width=0.05,color='b')
        fig.savefig(plotPath+filename+'_AggUB_AARanking.png')         
        return sorted_meanAAlis
    def GetAggregate_EntropyForAA_fromtxt (self,path,dictlist):
        """ This function takes dictionaries as text aggregates the dictionaries and returns the dict"""
        f = open(path+dictlist[0]+'.txt','r')
        my_dict0 = eval(f.read())
        f = open(path+dictlist[1]+'.txt','r')
        my_dict1 = eval(f.read())
        f = open(path+dictlist[2]+'.txt','r')
        my_dict2 = eval(f.read())
        f = open(path+dictlist[3]+'.txt','r')
        my_dict3 = eval(f.read())        
        
        aggdict={a:[] for a in my_dict3.keys()}
        for keys, value in my_dict0.iteritems():
            aggdict[keys].extend(value)
        for keys, value in my_dict1.iteritems():
            aggdict[keys].extend(value)           
        for keys, value in my_dict2.iteritems():
            aggdict[keys].extend(value)       
        for keys, value in my_dict3.iteritems():
            aggdict[keys].extend(value)                    
        return aggdict         
    def CalcSequenceEntropy(self,path):
        " This function is used to calculate sequence entropy using the fasta file and the homologous residues index list"
        f = FASTAnet.FASTAstructure(path, uniqueOnly=False, curate=False)
        aligIndex = f.findAlignedResidueIndices()
        seqs = f.getSequences()
        l=f.getSequenceLengths()
        count = 0
        myseq = []
        for s in seqs:
            newseq = []
            for q in s:         
                newseq.extend(q) 
            myseq.append(newseq)
        
        aalis=['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']    
        seqentropylist=[]
        for i in aligIndex:
            aadict={a:0 for a in aalis}
            totalfreq=0
            for s in myseq:
                key = s[i]
                aadict[key]=aadict[key]+1
                totalfreq=totalfreq+1  
            
            freqlis=[]   
            for keys,values in aadict.iteritems():
                if values > 0:
                    value=float(values)/totalfreq
                    freqlis.append(value*math.log(value,2))
            seqentropylist.append(-sum(freqlis))
        return seqentropylist
    def BoxPlotSeqEntropy(self,path,folderlist):
        r = ResEntropy()
        data0 = r.CalcSequenceEntropy(path+folderlist[0]+'/alignment.fasta')
        data1 = r.CalcSequenceEntropy(path+folderlist[1]+'/alignment.fasta')
        data2 = r.CalcSequenceEntropy(path+folderlist[2]+'/alignment.fasta')
        data3 = r.CalcSequenceEntropy(path+folderlist[3]+'/alignment.fasta')       
        
        data=[data0,data1,data2,data3]
        plt.boxplot(data)
        plt.boxplot(data,labels=folderlist)
        plt.ylim(-1,5)
        plt.show()
        fig = plt.figure(1, figsize=(20, 15))
        ax = fig.add_subplot(111)
        bp = ax.boxplot(data,labels=folderlist)
        fig.savefig(path+'Seqaggregate.png') 
    def HomoSeqEntropyGraph(self,path):
        """ This function reads a pdb file, gets the hoomolgous residues and plots entropy vs
            homologous sequence plot"""
        pdbRef=PDBnet.PDBstructure(path)
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        fasta = list(model.AsFASTA())
        homoindex=[r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        entropy=[r.GetAtoms()[0].tempFactor for r in residues if r.index in homoindex]
        seq=[s.name for s in residues if s.index in homoindex]
        #seq=[s for s,index in enumerate(fasta) if index in homoindex]
        for index,s in enumerate(fasta):
            if str(index) in homoindex:
                seq.append(s)
        x= range(len(seq))
        plt.title('Homologous Sequence vs Entropy Plot for Bound proteins')
        #plt.xticks(x, seq,rotation='vertical')
        plt.xlabel('Homologous Sequence')
        plt.ylabel('Entropy')
        plt.plot(entropy)
        plt.show()


    def PlotAA_Property_Entropy(self,entropydict, aminoacidSeq,title,xlabel,ylabel,filename,count,plotPath):
        """Input: entropydict: entropy value of each amino acid in the dictionary format
                   title, xlabel, ylabel for the plots
        Output: Plot for the increasing property of the amino acid There needs to be new method which is scattered plot to determine the relation"""
        xaxis = np.arange(0,20,1)
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle(title+filename,fontsize= 30)
        ax = fig.add_subplot(111)
        data = []
        for key in aminoacidSeq:
            data.extend([entropydict[key]])
        #data = entropydict.values()
        print data
        labels = aminoacidSeq
        #this is important so that the ticks are places at correct distance
        ax.set_xticks(xaxis)
        ax.set_xlabel(xlabel, fontsize= 20)
        ax.set_ylabel(ylabel, fontsize =20)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        ax.set_xticklabels(labels)
        ax.bar(xaxis,data,width=0.05,color='b')
        fig.savefig(plotPath+filename+'_Size.png') 
        print 'file has been made'
    def scatterPlotAA_Property(self,entropydict,csvfilePath,filename,count,outputPath):
        sizeDict = {}
        hydrophobicityDict = {}
        with open(csvfilePath,'rb') as csvfile:
            reader = csv.reader(csvfile,delimiter =',')
            for line in reader:
                hydrophobicityDict.update({line[0]:line[4]})
                sizeDict.update({line[0]:line[3]})
        #get the data for the plots
        xValues =[]
        ySize =[]
        yHydro = []
        for key, value in entropydict.iteritems():
            xValues.append(value)
            ySize.append(float(sizeDict[key]))
            yHydro.append(float(hydrophobicityDict[key]))
            
        mp = MyPlots()
        mp.simplescatterPlot(xValues, ySize, " Mean Entropy ", "Size of amino acids", "ScatterPlot between Mean Entropy and Size of "+filename,
                             outputPath+"Size"+filename,count)
        mp.simplescatterPlot(xValues, yHydro, " Mean Entropy ", "Hydrophobicity of amino acids", "Scatter Plot between Mean Entropy and Hydrophobicity of "+filename,
                            outputPath+"Hydrophobicity"+filename,count)        
    def SeqEntropyplot(self,path,activesitecoor):
        """ This function reads a pdb file,  plots entropy vs
            PDB sequence plot. This function also calls distligandresidue from the class active site
            to get the index of the active site"""
        pdb=PDBnet.PDBstructure(path)
        c = activeSite()
        activesiteindex = c.distligandresidue(path, (activesitecoor))
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        fasta = list(model.AsFASTA())
        entropy=[r.GetAtoms()[0].tempFactor for r in residues]
        seq=[s for s in fasta]
        
        #running average        
        weights = np.repeat(1.0, 10)/10
        sma = np.convolve(entropy, weights, 'valid')        
        #plot the entropy
        #plt.subplot(211)        
        #plt.plot([1,100,200,300])
        #plt.subplot(212)        
        x= range(len(seq))
        plt.title('Sequence vs Entropy Plot for Bound proteins')
        plt.ylim(-1,9)
        plt.xticks(x, seq,rotation='horizontal')
        plt.tick_params(labelsize=6)
        plt.xlabel('Homologous Sequence')
        plt.ylabel('Entropy')
        t=[1,100,200,300]
        plt.plot(sma,color="green",marker="o")
        plt.annotate('activesite', xy=(int(activesiteindex), sma[int(activesiteindex)]), xytext=(2,0.25),
                     arrowprops=dict(facecolor='black', shrink=0.05))
        plt.show()
        
        
    def SeqEntropyPlotfor2SequencesA(self,path1,path2,fastapath):
        
        """ This function reads a pdb file, gets plots entropy vs
            PDB sequence. This function also calls distligandresidue from the class active site
            to get the index of the active site. This function reads an alignment fasta file to get the alignment of bound and unbound
            proteins. It marks the entropy of the active site index and adds -0.01 to the black spaces in the alignment"""
        f = FASTAnet.FASTAstructure(fastapath, uniqueOnly=False, curate=False)
        seq =f.getSequences()
        alignSeq1 = []
        alignSeq2 = []
        for a in seq[1]:
            alignSeq1.extend(a)
        for a in seq[0]:
            alignSeq2.extend(a)
        index1 = [ i for i,x in enumerate(alignSeq1) if x == '-']
        index2 = [ i for i,x in enumerate(alignSeq2) if x == '-']
        pdb1=PDBnet.PDBstructure(path1)
        pdb2=PDBnet.PDBstructure(path2)
        c = activeSite()
        activesiteindex1 = c.distligandresidue(path1, (-23.716, -8.568, -8.136))
        m1 = pdb1.orderofmodels[0]
        model1=pdb1.GetModel(m1)
        residues1 = model1.GetResidues()
        seq1 = list(model1.AsFASTA())
        entropy1=[r.GetAtoms()[0].tempFactor for r in residues1]
        entropy1[activesiteindex1]='!'+str(entropy1[activesiteindex1])
        print  entropy1[activesiteindex1]
        print seq1[activesiteindex1]
        for i in index1:
            entropy1.insert(i,-0.01)
        for i,j in enumerate(entropy1):
            if '!' in str(j):
                entropy1[i] = float(str(entropy1[i].replace('!','')))
                asIndex1 = i
        print  entropy1[activesiteindex1]
        print alignSeq1[activesiteindex1]
        activesiteindex2 = c.distligandresidue(path2, (115.429,56.667,48.841))
        m2 = pdb2.orderofmodels[0]
        model2=pdb2.GetModel(m2)
        residues2 = model2.GetResidues()
        seq2 = list(model2.AsFASTA())
        entropy2=[r.GetAtoms()[0].tempFactor for r in residues2]  
        entropy2[activesiteindex2]='!'+str(entropy1[activesiteindex2])
        for i in index2:
            entropy2.insert(i,-0.01)        
        for i,j in enumerate(entropy2):
            if '!' in str(j):
                entropy2[i]=float(str(entropy2[i].replace('!','')))
                asIndex2 = i        
        #difference =[abs(s1-s2) for s1,s2 in zip(entropy1,entropy2)]
        #running average        
        weights = np.repeat(1.0, 10)/10
        sma1 = np.convolve(entropy1, weights, 'valid')
        sma2 = np.convolve(entropy2, weights, 'valid')
        #plot the entropy
        x= range(len(alignSeq1))
        plt.subplot(111)
        plt.title('Sequence vs Running Average Entropy Plot for Bound Unbound GlucarateDehydratase')
        plt.ylim(-1,6)
        plt.xticks(x, alignSeq2,rotation='horizontal')
        plt.tick_params(labelsize=10)
        plt.xlabel('Homologous Sequence')
        plt.ylabel('Entropy')
        plt.plot(sma1,'g-')
        plt.annotate('activesiteB', xy=(asIndex1, sma1[asIndex1]), xytext=(asIndex1,sma1[asIndex1]+2),
                     arrowprops=dict(facecolor='black', shrink=0.05))        
                
        #plt.plot([1,100,200,300])
        #plt.subplot(212)        
        #x= range(len(seq2))
        #plt.title('Sequence vs Running Average Entropy Plot for UnBound proteins')
        #plt.ylim(-1,3)
        #plt.xticks(x, seq2,rotation='horizontal')
        #plt.tick_params(labelsize=6)
        #plt.xlabel('Homologous Sequence')
        #plt.ylabel('Entropy')
        plt.plot(sma2,'b-')
        #plt.plot(difference,'r-')
        plt.annotate('activesiteUB', xy=(asIndex2, sma2[asIndex2]), xytext=(asIndex2,sma2[asIndex2]+1),
                     arrowprops=dict(facecolor='red', shrink=0.05))
        plt.show() 
    def SeqEntropyPlotfor2SequencesB(self,path1,path2,fastapath,title,activeSitecoorB,activeSitecoorUB):
        #this function is the second version of SeqEntropyPlotfor2SequencesA because that thing was clearly not working
            """ This function reads a pdb file, gets plots entropy vs
                PDB sequence. This function also calls distligandresidue from the class active site
                to get the index of the active site. This function reads an alignment fasta file to get the alignment of bound and unbound
                proteins. It marks the entropy of the active site index and adds -0.01 to the black spaces in the alignment.
                This fucntion takes 2 sequences as input and gets the entropy of only aligned amino acids. It plots the aligned sequence
                and the entropy"""                   
            f = FASTAnet.FASTAstructure(fastapath, uniqueOnly=False, curate=False)
            seq =f.getSequences()
            ungappedindex = f.getStrictlyUngappedPositions()
            alignSeq1 = []
            alignSeq2 = []
            if len(seq)<2:
                print "The fasta file not found"
                
            for a in seq[1]:
                alignSeq1.extend(a)
            for a in seq[0]:
                alignSeq2.extend(a) 
            pdb1=PDBnet.PDBstructure(path1)
            pdb2=PDBnet.PDBstructure(path2)
            modelname = pdb1.GetModelNames()         
            c = activeSite()          
            activesiteindex1 = c.distligandresidue(path1, activeSitecoorB)
            m1 = pdb1.orderofmodels[0]
            model1=pdb1.GetModel(m1)
            residues1 = model1.GetResidues()
            seq1 = list(model1.AsFASTA())
            entropy1=[r.GetAtoms()[0].tempFactor for r in residues1]
            entropy=[r.GetAtoms()[0].tempFactor for r in residues1]
            entropy1[activesiteindex1]='!'+str(entropy1[activesiteindex1])
            print seq1[activesiteindex1]
            print activesiteindex1
            ungappedentropy1=[]
            offset = 1
            for i in enumerate(alignSeq1): 
                if i[1] =='-':
                    entropy1.insert(i[0]+offset,-1)
                    offset = offset+1
            for i in ungappedindex: 
                ungappedentropy1.append(entropy1[i])
            for i,j in enumerate(ungappedentropy1):
                if '!' in str(j):
                    ungappedentropy1[i] =  float(str(ungappedentropy1[i].replace('!','')))
                    asIndex1 = i
                    
                    
                    
            activesiteindex2 = c.distligandresidue(path2, activeSitecoorUB)            
            m2 = pdb2.orderofmodels[0]
            model2=pdb2.GetModel(m2)
            residues2 = model2.GetResidues()
            seq2 = list(model2.AsFASTA())
            entropy2=[r.GetAtoms()[0].tempFactor for r in residues2]  
            entropy2[activesiteindex2]='!'+str(entropy2[activesiteindex2])
            print seq2[activesiteindex2]
            print activesiteindex2
            ungappedentropy2=[]
            offset = 1
            for i in enumerate(alignSeq2): 
                if i[1] =='-':
                    entropy2.insert(i[0]+offset,-0.01)
                    offset = offset+1        
            for i in ungappedindex: 
                ungappedentropy2.append(entropy2[i])              
            for i,j in enumerate(ungappedentropy2):
                if '!' in str(j):
                    ungappedentropy2[i]=float(str(ungappedentropy2[i].replace('!','')))
                    asIndex2 = i         
            weights = np.repeat(1.0, 10)/10
            sma1 = np.convolve(ungappedentropy1, weights, 'valid')
            sma2 = np.convolve(ungappedentropy2, weights, 'valid')
            #difference =[max(s1,s2)-min(s1,s2) for s1,s2 in zip(ungappedentropy1,ungappedentropy2)] 
            difference = [s1-s2 for s1, s2 in zip(ungappedentropy1,ungappedentropy2)]
            #plot the entropy
            x= range(len(alignSeq1))
            plt.subplot(111)
            plt.title('Sequence vs Running Average Entropy Plot for Bound Unbound '+title)
            plt.ylim(-1,10)
            xticks = [alignSeq1[i] for i in ungappedindex]
            plt.xticks(x, xticks,rotation='horizontal')
            plt.tick_params(labelsize=8)
            plt.xlabel('Homologous Sequence',fontsize = 20)
            plt.ylabel('Entropy',fontsize = 20)
            #plt.plot(ungappedentropy1,'b-')
            plt.annotate('activesiteB', xy=(asIndex1, sma1[asIndex1]), xytext=(asIndex1,sma1[asIndex1]+3),
            
            
            
            
                         arrowprops=dict(facecolor='black', shrink=0.05))
            #plt.plot(ungappedentropy2,'g-')
            plt.plot(difference,'-r',label='Difference')
            plt.annotate('activesiteUB', xy=(asIndex2, sma2[asIndex2]), xytext=(asIndex2,sma2[asIndex2]+3),
                         arrowprops=dict(facecolor='red', shrink=0.05))
            #line_up = plt.plot([1,2,3], label='Bound')
            #line_down = plt.plot([3,2,1], label='UnBound')
            plt.legend()
            plt.show()
            q = pdb1._FastaPdbMatch(fastapath)
            pdb1.Map2Protein(path1+'3.pdb',difference,0,fastapath) 
            
    def PlotContactEntropyvsSequenceEntropy(self, contactEntropy, SequenceEntropy,homosequence,title):
        x= range(len(homosequence))
        plt.subplot(111)
        plt.title(title)
        plt.ylim(0,10)
        plt.xticks(x, homosequence,rotation='horizontal')
        plt.tick_params(labelsize=8)
        plt.xlabel('Homologous Sequence',fontsize = 20)
        plt.ylabel('Entropy Values',fontsize = 20)
        plt.plot(contactEntropy,'b-')
        plt.plot(SequenceEntropy,'g-')
        line_up = plt.plot([1,2,3], label='Contact Entropy')
        line_down = plt.plot([3,2,1], label='Sequence Entropy')
        plt.legend()
        plt.show()        
        
class activeSite:
    def __init__(self, pdbfile=None):

        if pdbfile == None:
            self.pdb = None
        else:
            self.pdb = pdbfile
    
        self.cent   = None 
        
    def distanceActiveSiteIndices(self,pdb,asindex,compPath):
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        
        homoindex=[r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        entropy=[r.GetAtoms()[0].tempFactor for r in residues if r.index in homoindex]
        ascent=[a.Centroid() for a in residues if a.index in asindex]
        distance=[]
        #get the centroid of all the residues
        for res in residues:
            if res.index in homoindex:
                centroid = res.Centroid()
                dis=min([centroid.DistanceTo(a) for a in ascent])
            distance.append(min([centroid.DistanceTo(a) for a in ascent]))
        tup=[[d,e] for d,e in zip(distance,entropy)]
        np.save(compPath, tup)
        plt.title('ActiveSiteDistancePlot')
        plt.xlabel('DistanceFromActiveSite')
        plt.ylabel('Entropy')
        plt.plot(distance,entropy,"o")
        #plt.show()
        print 'done!'
    
    
    def distFromLigands_entropy(self,pdb,compPath,title,rand,ligandcoor):
        """This function assumes that mg and mn are the active sites of the protein. This function access a new function in PDBnet
               def DistanceToligand(self, ligandcord):
                   return distance.dist(self.GetPosition(),ligandcord)."""
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        homoindex=[r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        entropy=[r.GetAtoms()[0].tempFactor for r in residues if r.index in homoindex]
        enolaseB=(26.634 ,74.962 ,12.870)
        gdB=(-23.716, -8.568, -8.136)
        mrB =(-20.292,60.052 ,55.442)
        mleB =(47.778,31.504,34.784)
        enolaseUB=(20.212,17.290,27.123)
        gdUB=(-17.144,-17.786,-3.096)
        mrUB=(-19.509,60.979,56.437)
        mleUB=(-7.516,16.764,13.258)
        distance=[]
        activeSiteindex = activeSite.distligandresidue(self, compPath+'.pdb',ligandcoor)
        coordinates = [res.Centroid().GetPosition() for res in residues]
        randomsite = random.choice(coordinates)
        activesitecoor = coordinates[int(activeSiteindex)]
        #get the centroid of all the residues
        for res in residues:
            if res.index in homoindex:
                centroid = res.Centroid()
                dis = centroid.DistanceToligand(activesitecoor)
                distance.append(dis)
        Proteintuple=[[d,e] for d,e in zip(distance,entropy)]
        tup=[[d,e] for d,e in zip(distance,entropy) if d<10]
        randomtup = []
        x = [ a[0] for a in tup]
        y = [ a[1] for a in tup]        

        #np.save(compPath, tup)
        plt.title('Active Site for'+title)
        plt.ylim(0,8)
        plt.xlabel('DistanceFromActiveSite')
        plt.ylabel('Entropy')

        
        for i in range(0,len(x)):
            randomtup.append(random.choice(Proteintuple))
        rx = [ a[0] for a in randomtup]
        ry = [ a[1] for a in randomtup]
        if rand == True:
            f = open(compPath+'random2.csv', 'w')
            writer = csv.writer(f,delimiter='\t') 
            for t in randomtup:
                writer.writerow(t)
            f.close()            
            plt.plot(rx,ry,"o")
            plt.show()
        else:      
            f = open(compPath+'.csv', 'w')
            writer = csv.writer(f,delimiter='\t') 
            for t in tup:
                writer.writerow(t)
            f.close()            
            plt.plot(x,y,"o")
            plt.show()        
        
    def distligandresidue(self, path,ligandcoordinates): 
        """This function takes a PDB file and ligand coordinates as input and finds the index of the residue 
        which is closest to the ligand"""
        pdb=PDBnet.PDBstructure(path)
        if pdb.orderofmodels ==[]:
            ch = pdb.GetChainNames()
            chain = pdb.GetChain(ch[0])
            residues = chain.GetResidues()
        else:
            m = pdb.orderofmodels[0]
            model=pdb.GetModel(m)  
            residues = model.GetResidues()
        
        residuedict ={x.index:0 for x in residues}
        distance = []
        for r in residues:
                centroid = r.Centroid()
                dis = centroid.DistanceToligand(ligandcoordinates)
                residuedict[r.index] = dis
                distance.append(dis)
        distance.sort()
        for k,v in residuedict.iteritems():
            if v==distance[0]:
                activesiteindex = k
        return int(activesiteindex)
    def DrawAggScatterPlot(self,path,dictlist):
        tup0=np.load(path+dictlist[0]+'.npy')
        tup1=np.load(path+dictlist[1]+'.npy')
        tup2=np.load(path+dictlist[2]+'.npy')
        tup3=np.load(path+dictlist[3]+'.npy')
        distance=[]
        entropy=[]
        for a,b in tup0:
            distance.append(a)
            entropy.append(b)
        for a,b in tup1:
                distance.append(a)
                entropy.append(b)       
        for a,b in tup2:
                distance.append(a)
                entropy.append(b)        
        for a,b in tup3:
                distance.append(a)
                entropy.append(b)
        plt.title('AggregateActiveSiteDistancePlot')
        plt.xlabel('DistanceFromActiveSite')
        plt.ylabel('Entropy')
        plt.plot(distance,entropy,"o")
        plt.show()
        
class PDBChainComp:
    def ChainCompfor2(self,path):
        """ Compare if the two chains in a PDB are equal if they are Calculate the FCM matrix to see the structural difference.
            This function """
        output = open('/home/qazi/mydata/LabWork/Documents/SimilarChains2.txt','w')
        files = os.listdir(path)
        for f in files:
            print f
            p = PDBnet.PDBstructure(path+f)
            chainNames = p.GetChainNames()
            seqdict={x:[] for x in chainNames}
            for c in chainNames:
                ch = p.GetChain(c)
                seq=ch.AsFASTA()
                seqdict[c].append(seq)
            if len(seqdict)>1:
                
                for a, b in itertools.combinations(seqdict, 2):
                    flag=0
                    count =0
                    value = seqdict[a] == seqdict[b] 
                    print a,b,value
                    if value==True:
                        m1 = p.ContactMatrix(a)
                        m2 =p.ContactMatrix(b)
                        aggmtrix=np.zeros((len(m1), len(m1[0])))
                        for i in range(len(m1)):
                            for j in range(len(m1[0])):
                                aggmtrix[i][j] = (m1[i][j] + m2[i][j])/2
                        for element in aggmtrix.flat:
                            if 0.0 < element < 1.0:
                                flag=1
                                count =count+1
                                #break
                        print count
                        if flag==1:
                            output.write(f+" has chain "+a+" and chain "+b+" with different structures but same sequence with "+str(count)+"\n")
                            output.flush()
                            flag=0
            else:
                output.write(f+" has only one chain\n")
                output.flush()
        
        output.close()
        print 'done!'
                
        
    def ChainComp(self,path):
        
        
        output = open('/home/qazi/mydata/LabWork/Documents/SimilarChains4.txt','w')
        files = os.listdir(path)
        for f in files:
            print f            
            p = PDBnet.PDBstructure(path+f)
            chainNames = p.GetChainNames()
            seqdict={x:[] for x in chainNames}
            for c in chainNames:
                ch = p.GetChain(c)
                seq=ch.AsFASTA()
                seqdict[c]=seq
            newdict = {}
            if len(seqdict)>1:
                for key, value in seqdict.items():
                    newdict.setdefault(value,set()).add(key)
                keys=[values for key, values in newdict.items() if len(values) > 1]
                print keys
                
                for keyset in keys:
                    count = 0
                    countones = 0 
                    keylist=list(keyset)
                    aggmatrix=p.ContactMatrix(keylist[0])
                    for index in range(1,len(keylist)):
                        m = p.ContactMatrix(keylist[index])
                        aggmatrix=np.add(aggmatrix,m)
                    fcm=np.divide(aggmatrix,len(keyset))
                    plt.matshow(fcm)
                    plt.show()
                    #np.savetxt(path+f+'FCM', fcm)
                    for element in fcm.flat:
                        if 0.0 < element < 1.0:
                            flag=1
                            count =count+1
                            #break
                        if abs(element-1.0)==0:
                            countones =countones+1
                    print count
                    if flag==1:
                        filename =f+' '
                        output.write(filename+" ".join([k for k in keylist])+" "+str(count)+" "+str(countones)+" \n")
                        output.flush()
                        flag=0                        
                    
            else:
                output.write(f+" has only one chain\n")
                output.flush()
        
        output.close()
        print 'done!' 
    def Entropyfor1PDB(self,path):
        p = PDBnet.PDBstructure(path)
        chainNames = p.GetChainNames()
        entropy=[]
        seqdict={x:[] for x in chainNames}
        for c in chainNames:
            ch = p.GetChain(c)
            seq=ch.AsFASTA()
            seqdict[c]=seq
        newdict = {}
        if len(seqdict)>1:
            
            for key, value in seqdict.items():
                newdict.setdefault(value,set()).add(key)
            keys=[values for key, values in newdict.items() if len(values) > 1]
            print keys
            for keyset in keys:
                count = 0
                countones = 0 
                keylist=list(keyset)
                aggmatrix=p.ContactMatrix(keylist[0])
                for index in range(1,len(keylist)):
                    m = p.ContactMatrix(keylist[index])
                    aggmatrix=np.add(aggmatrix,m)
                fcm=np.divide(aggmatrix,len(keyset))
                e=ResEntropy()
                entropy = e.calcEntropy(fcm)
        return entropy
    def densityPlotforPDB(self,path,pathPDB):
        
        files=os.listdir(path)
        lfiles=[f.lower() for f in files]
        pdbs=os.listdir(pathPDB)
        lpdbs=[f.lower() for f in pdbs]
        n=[]
        
        for f in lfiles:
            n.append(f[:4])
        s=set(n)
        subnames=list(s)
        entropy=[]
        for name in subnames:
            for i in lpdbs:
                if i[:4]==name:
                    print i
                    newlist = self.Entropyfor1PDB(pathPDB+i)
                    entropy=entropy+newlist
        density = gaussian_kde(entropy)
        xs = np.linspace(min(entropy),max(entropy))
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        plt.plot(xs,density(xs))
        plt.show()                    
        print names
    def comparechainsforsimilarity(self,folderpath1,folderpath2):
        """This function takes chains from one folder and compares it to the chains(PDB)from another folder to get ones with similar chains"""
        boundpdb = os.listdir(folderpath1)
        unboundpdb = os.listdir(folderpath2)
        
        for bound in boundpdb:
            pdb1=PDBnet.PDBstructure(folderpath1+bound)
            chn1=bound[4].upper()
            chain1= pdb1.GetChain(chn1)
            seq1 = chain1.AsFASTA()
            for unbound in unboundpdb:
                counter = 0
                pdb2=PDBnet.PDBstructure(folderpath2+unbound)
                chn2=unbound[4].upper()
                chain2= pdb2.GetChain(chn2)
                seq2 = chain2.AsFASTA()
                for s1,s2 in zip(seq1,seq2):
                    if s1 == s2:
                        counter = counter+1
                print bound, ': ',len(seq1), unbound,': ',len(seq2),' have ',counter         
                #if seq1==seq2:
                    #print bound, 'and', unbound
class simulatedata:
    def __init__(self):
        self.pdbs = None
    def SimplePlot(self,yValues,xlables,outputPath):
        xaxis = np.arange(1,6,1)
        fig = plt.figure(1, figsize=(30, 20))
        fig.suptitle('Randomized Values in FCM for Bound',fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Families', fontsize= 20)
        ax.set_ylabel('Number of Randomized values in FCM', fontsize =20)
        ax.set_xticks(xaxis)
        ax.set_xticklabels(xlables,rotation = 'vertical')
        ax.tick_params(axis = 'x', labelsize=20)
        ax.tick_params(axis = 'y', labelsize=30)
        pt = ax.bar(xaxis,yValues,width=0.05,color='r')      
        fig.savefig(outputPath+'.png')        
        
    def simulatenormaldata(self,pdbpath,fcm,ligandcoor,csvpath):
        """Input: path of the pdb file that needs to be used to get the active site
                  fcm : fcm[0] is the frequencey contact matrix and fcm[1] is the fastaresidue homologs
                  ligandcoordinates: are the coordinates of the ligand that are used to find the active site
           Output: a csv file with FCm values of the residues that are at a distance of 10 A from the active site 
                   and other normal distribution with same mean and different standard deviation truncated between 0 and 1
        """
        acsite = activeSite()
        activeSiteindex = acsite.distligandresidue(pdbpath,ligandcoor)
        distance={}
        pdb = PDBnet.PDBstructure(pdbpath)
        
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        coordinates = [res.Centroid().GetPosition() for res in residues]
        activesitecoor = coordinates[int(activeSiteindex)]
        #tempresidues = []
        #for i, pos in enumerate(fcm[1]):
            #res = model.GetResidueByPosition(pos)
            #tempresidues.append(res.index)         
        homoindex=[ r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        distance = {x:0 for x in homoindex}
        for res in residues:
            if res.index in homoindex:
                centroid = res.Centroid()
                dis = centroid.DistanceToligand(activesitecoor)
                distance[res.index]=dis
        
        resindex=[int(index) for index,d in distance.iteritems() if d<10]
        #print resindex
        matrixindex = [homoindex.index(str(i)) for i in resindex]
        combos = [combo for combo in itertools.combinations(matrixindex, 2)]
        f = file(csvpath,'wb')
        writer = csv.writer(f, delimiter = '\t')
        data =[] 
        for comb in combos:
            data.append(fcm[0][comb[0]][comb[1]])
        filterdata = filter(lambda a :a != 0.0, data) 
        #print filterdata
        mycsv = []
        mycsv.append(filterdata)
        shape = len(filterdata)
        mean = np.mean(filterdata)
        lower, upper = 0.0, 1.0
        sigmavalues = np.arange(0,1,0.05)
        for sigma in sigmavalues:
            #normaldataobj = stats.truncnorm((lower - mean) / sigma, (upper - mean) / sigma, loc=mean, scale=sigma)
            #normaldata = normaldataobj.rvs(shape)
            normaldata =[]
            for i in filterdata:
                myvalue = min(max(gauss(float(i),float(sigma)),0.0),1.0)
                normaldata.append(myvalue)
            #print normaldata
            mycsv.append(normaldata)           
            #plt.plot(normaldata.rvs(shape),'o')
            #plt.show()
        mycsv2 = map(list, zip(*mycsv))
        for row in mycsv2:
            writer.writerow(row)
        f.close()
        print 'finish'
    def simulate100normaldata(self,pdbpath,fcm,ligandcoor,csvpath):
            """Input: path of the pdb file that needs to be used to get the active site
                      fcm : fcm[0] is the frequencey contact matrix and fcm[1] is the fastaresidue homologs
                      ligandcoordinates: are the coordinates of the ligand that are used to find the active site
               Output: a csv file with FCM values of the residues that are at a distance of 10 A from the active site 
                       and other 100 normal distribution with same mean and different standard deviation truncated between 0 and 1
            """
            acsite = activeSite()
            activeSiteindex = acsite.distligandresidue(pdbpath,ligandcoor)
            distance={}
            pdb = PDBnet.PDBstructure(pdbpath)
            if pdb.orderofmodels ==[]:
                ch = pdb.GetChainNames()
                chain = pdb.GetChain(ch[0])
                residues = chain.GetResidues()
            else:
                m = pdb.orderofmodels[0]
                model=pdb.GetModel(m)  
                residues = model.GetResidues()
            coordinates = [res.Centroid().GetPosition() for res in residues]
            activesitecoor = coordinates[int(activeSiteindex)]
            #homoindex = [ r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
            homoindex = [str(v) for v in fcm[1]]
            distance = {x:0 for x in homoindex}
            for res in residues:
                if res.index in homoindex:
                    centroid = res.Centroid()
                    dis = centroid.DistanceToligand(activesitecoor)
                    distance[res.index]=dis
            
            resindex=[int(index) for index,d in distance.iteritems() if d<10]
            matrixindex = [homoindex.index(str(i)) for i in resindex]
            #matrixindex = [homoindex.index(i) for i in resindex]
            combos = [combo for combo in itertools.combinations(matrixindex, 2)]
            data =[] 
            for comb in combos:
                data.append(fcm[0][comb[0]][comb[1]])
            filterdata = filter(lambda a :a != 0.0, data)
            f = file(csvpath+'.csv','wb')
            writer = csv.writer(f, delimiter = '\t')
            for row in filterdata:
                writer.writerow([row])
            f.close()
            shape = len(filterdata)
            sigmavalues = np.arange(0,1,0.05)
            count = 0          
            for sigma in sigmavalues:
                print 'sigma ',sigma
                mycsv=[]
                for number in range(0,100):
                    random.seed() 
                    normaldata =[]
                    for i in filterdata:
                        myvalue = min(max(gauss(float(i),float(sigma)),0.0),1.0)
                        #print i ,' ',myvalue
                        normaldata.append(myvalue)
                    mycsv.append(normaldata)
                #plt.plot(filterdata,color ='r')
                #plt.plot(normaldata,color ='b')
                #plt.show()
                #plot the two distributions
                #count = count+1
                #fig = plt.figure(count, figsize=(30, 20))
                #ax = fig.add_subplot(111)
                #ax.tick_params(axis = 'x', labelsize=30)
                #ax.tick_params(axis = 'y', labelsize=30)
                #ax.plot(filterdata,color ='r')
                #ax.plot(normaldata,color ='b')
                #fig.savefig(csvpath+'_'+str(sigma)+'.png')
                
                
                
                f2 = file(csvpath+'_'+str(sigma)+'.csv','wb')
                writer1 = csv.writer(f2, delimiter = '\t')                    
                transposecsv = map(list, zip(*mycsv))
                for row in transposecsv:
                    writer1.writerow(row)
                f2.close()
                print f2
    def plotPvalue(self,inpath,outpath,outnames):
        
        xaxis = np.arange(0,1,0.05)
        count = 0
        with open(inpath,'rb') as  csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            for row in reader:
                #plt.title('P Value vs Sigma '+outnames[count])
                #plt.xlabel('Sigma')
                #plt.ylabel('P value')                
                #plt.plot(xaxis,row)
                #plt.show()
                fig = plt.figure(count, figsize=(30, 20))
                fig.suptitle('P Value vs Sigma '+outnames[count],fontsize= 30)
                ax = fig.add_subplot(111)
                ax.set_xlabel('Sigma', fontsize= 20)
                ax.set_ylabel('P value', fontsize =20)
                ax.tick_params(axis = 'x', labelsize=30)
                ax.tick_params(axis = 'y', labelsize=30)
                pt = ax.plot(xaxis,row)
                fig.savefig(outpath+'/'+str(count)+'.png')
                count = count+1
                print 'this'
    def plotPvaluesfromIndividualcsv(self,inpath,outpath):
        filenames = os.listdir(inpath)
        xaxis = np.arange(0,1,0.05)
        count = 0
        for f in filenames: 
            print f[:-4]
            with open(inpath+f,'rb') as  csvfile:
                reader = csv.reader(csvfile, delimiter = ',')
                values=[]
                newvalues=[]
                for row in reader:                    
                    row2 = list(map(float,row))
                    # this is to get the bar plot for the values that are less than 0.05
                    lessthanval = sum(i<0.05 for i in row2)
                    newvalues.append(lessthanval)
                    #values.append(row2)
                    #filterrow = filter(lambda a:a<0.05,row2)
                    #values.append(filterrow)
                #nonempty =[]
                #for x in values:
                    #if not x:
                        #nonempty.append([0]*100)
                    #else:
                        #nonempty.append(x)
                #transposedvalues= zip(*nonempty)
                transposedvalues = map(list,zip(*values))
                fig = plt.figure(count, figsize=(30, 20))
                fig.suptitle('P Value vs Sigma '+f,fontsize= 30)
                ax = fig.add_subplot(111)
                ax.set_xlabel('Sigma', fontsize= 20)
                ax.set_ylabel('P value', fontsize =20)
                ax.tick_params(axis = 'x', labelsize=30)
                ax.tick_params(axis = 'y', labelsize=30)
                
                #ax.boxplot(values,labels=xaxis)
                #you just need the numbre of values that are less than 0.05 
                ax.bar(xaxis,newvalues,width=0.01,color='r')
                #for y in transposedvalues:
                    #ax.set_ylim([0.0,0.05])
                    #ax.bar(xaxis,y,width=0.01,color='r')
                    #ax.plot(xaxis,y,'o')
                    
                    
                fig.savefig(outpath+f[:-4]+'.png')               
                count = count+1
    def plot_sum_Pvalues_lessthan(self,inpath,outpath):
        """this is a cleaner version of plotPvaluesfromIndividualcsv. 
        It takes the csv for every fileas input and then calculates the number of values that 
        are less than 0.05 for each sigma
        """
        filenames = os.listdir(inpath)
        xaxis = np.arange(0,1,0.05)
        yaxis = np.arange(0,100,1)
        count = 0
        for f in filenames: 
            print f[:-4]
            with open(inpath+f,'rb') as  csvfile:
                reader = csv.reader(csvfile, delimiter = ',')
                newvalues=[]
                for row in reader:                    
                    row2 = list(map(float,row))
                    lessthanval = sum(i<0.05 for i in row2)
                    print lessthanval
                    newvalues.append(lessthanval)
                fig = plt.figure(count, figsize=(30, 20))
                fig.suptitle('Number of P values less than 0.05 vs Sigma '+f,fontsize= 30)
                ax = fig.add_subplot(111)
                ax.set_xlabel('Sigma', fontsize= 20)
                ax.set_ylabel('Number of P values less than 0.05', fontsize =20)
                ax.tick_params(axis = 'x', labelsize=30)
                ax.tick_params(axis = 'y', labelsize=30)
                ax.set_ylim([1,100])
                ax.bar(xaxis,newvalues,width=0.01,color='r')        
                fig.savefig(outpath+f[:-4]+'.png')               
                count = count+1
    def plot_sum_Pvalues_lessthan_for_2_datasets(self,inpath,outpath):
        """This function is same as that of plot_sum_Pvalues_lessthan but reads two csv files at once instead of reading things from the folder 
        one by one
        """
        filenames = os.listdir(inpath)
        xaxis = np.arange(0,1,0.05)
        yaxis = np.arange(0,100,1)
        count = 0
        print filenames[0]
        with open(inpath+filenames[0],'rb') as  csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            newvalues=[]
            for row in reader:                    
                row2 = list(map(float,row))
                lessthanval = sum(i<0.05 for i in row2)
                print lessthanval
                newvalues.append(lessthanval)
        print filenames[1]
        with open(inpath+filenames[1],'rb') as  csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            newvalues2=[]
            for row in reader:                    
                row2 = list(map(float,row))
                lessthanval2 = sum(i<0.05 for i in row2)
                print lessthanval2
                newvalues2.append(lessthanval2)        
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle('Number of P values less than 0.05 vs Sigma combined',fontsize= 30)
        ax = fig.add_subplot(111)
        w = 0.3
        ax.set_xlabel('Sigma', fontsize= 20)
        ax.set_ylabel('Number of P values less than 0.05', fontsize =20)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        ax.set_ylim([1,100])
        rect1 = ax.bar(xaxis-0.01,newvalues,width=0.01,color='r')
        rect2 = ax.bar(xaxis,newvalues2,width=0.01,color='b')
        ax.legend((rect1[0],rect2[0]),('Unbound','Bound'))
        fig.savefig(outpath+'combined.png')               
        count = count+1     
    def visualizeMatrix(self, matrix,outputpath):
        #plt.imshow(matrix,interpolation='nearest',cmap='hot',vmin = 0, vmax = 1)
        #plt.show()
        fig = plt.figure(0, figsize=(30, 20))
        fig.suptitle('Matrix',fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Amino Acids', fontsize= 20)
        ax.set_ylabel('Amino Acids', fontsize =20)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        pt = ax.imshow(matrix,interpolation='nearest',cmap='hot',vmin = 0, vmax = 1)        
        fig.savefig(outputpath)
        
    def create_SD_FCM (self,pdbpath,fcm,ligandcoor,csvpath,matrixpath):
        """This function siumates the data 100 times for every sigma. It also calls the function visualize matrix to create an image of a matrix 
        """
        acsite = activeSite()
        activeSiteindex = acsite.distligandresidue(pdbpath,ligandcoor)
        s = simulatedata()
        distance={}
        pdb = PDBnet.PDBstructure(pdbpath)
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        coordinates = [res.Centroid().GetPosition() for res in residues]
        activesitecoor = coordinates[int(activeSiteindex)]       
        homoindex=[ r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        distance = {x:0 for x in homoindex}
        for res in residues:
            if res.index in homoindex:
                centroid = res.Centroid()
                dis = centroid.DistanceToligand(activesitecoor)
                distance[res.index]=dis
        
        resindex=[int(index) for index,d in distance.iteritems() if d<10]
        matrixindex = [homoindex.index(str(i)) for i in resindex]
        combos = [combo for combo in itertools.permutations(matrixindex, 2)]
        totalCombos = [combo for combo in itertools.permutations(range(len(fcm[0])), 2)]
        #create SDFCM
        data =[] 
        for comb in combos:
            data.append(fcm[0][comb[0]][comb[1]])
        filterdata = filter(lambda a :a != 0.0, data)
        f = file(csvpath+'.csv','wb')
        writer = csv.writer(f, delimiter = '\t')
        for row in filterdata:
            writer.writerow([row])
        f.close()
        shape = len(filterdata)
        mean = np.mean(filterdata)
        lower, upper = 0.0, 1.0
        sigmavalues = np.arange(0,1,0.05)
        for sigma in sigmavalues:
            mycsv=[]
            for number in range(0,100):
                random.seed() 
                normaldata =[]
                for i in filterdata:
                    myvalue = min(max(gauss(float(i),float(sigma)),0.0),1.0)
                    normaldata.append(myvalue)
                mycsv.append(normaldata)                
            temp = 0
            mymatrix = fcm[0]
            for c in combos:
                if mymatrix[c[0]][c[1]] != 0.0:
                    mymatrix[c[0]][c[1]] = normaldata[temp]
                    temp = temp+1
            for x in totalCombos:
                if x not in combos:
                    mymatrix[x[0]][x[1]] = 0.0
            s.visualizeMatrix(mymatrix, matrixpath+str(sigma)+'.png')
            f2 = file(csvpath+'_'+str(sigma)+'.csv','wb')
            writer1 = csv.writer(f2, delimiter = '\t')                    
            transposecsv = map(list, zip(*mycsv))
            for row in transposecsv:
                writer1.writerow(row)
            f2.close()
    def FCM_Values_ClosetoActiveSite(self,pdbpath,fcm,ligandcoor):
        """This function counts the number of values in the FCM (also the values) that are close to the active site and are not equal to 0.0
        and plot the values for all the sub families
        """
        acsite = activeSite()
        activeSiteindex = acsite.distligandresidue(pdbpath,ligandcoor)
        distance={}
        pdb = PDBnet.PDBstructure(pdbpath)
        if pdb.orderofmodels ==[]:
            ch = pdb.GetChainNames()
            chain = pdb.GetChain(ch[0])
            residues = chain.GetResidues()
        else:
            m = pdb.orderofmodels[0]
            model=pdb.GetModel(m)  
            residues = model.GetResidues()
        coordinates = [res.Centroid().GetPosition() for res in residues]
        activesitecoor = coordinates[int(activeSiteindex)]       
        #homoindex=[ r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        homoindex = [str(v+1) for v in fcm[1]]
        distance = {x:0 for x in homoindex}
        for res in residues:
            if res.index in homoindex:
                centroid = res.Centroid()
                dis = centroid.DistanceToligand(activesitecoor)
                distance[res.index]=dis
        
        resindex=[int(index) for index,d in distance.iteritems() if d<10]
        matrixindex = [homoindex.index(str(i)) for i in resindex]
        combos = [combo for combo in itertools.combinations(matrixindex, 2)]
        data =[] 
        for comb in combos:
            data.append(fcm[0][comb[0]][comb[1]])
        filterdata = filter(lambda a :a != 0.0, data)
        return filterdata
        #return (len(filterdata))
        
    def PlotFCM_2_Distribution(self, data1,data2,outputPath,count):
        """This function gets two data distributions and plots the distributions and their difference"""
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle('Distribution of FCM Bound and Unbound',fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_ylabel('FCM values', fontsize= 20)       
        ax = fig.add_subplot(111)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        dist1 = ax.plot(data1,color ='r')
        dist2 = ax.plot(data2,color ='b')
        ax.legend((dist1[0],dist2[0]),('Unbound','Bound'))
        fig.savefig(outputPath+'.png')        
        
         
        
    def BarPlots_FCMDistribution(self,data,outputPath,count,plotname):
        """This function generates histograms for the number of times a value appears in the FCM"""
        fig = plt.figure(count, figsize=(30, 20))
        xaxis = np.arange(0,1,0.05)
        fig.suptitle('Bar chart for FCM values close to Active site for '+plotname,fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel('FCM values', fontsize= 40)    
        ax.set_ylabel('Number of FCM values of a particular value', fontsize= 40)    
        ax = fig.add_subplot(111)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        dist1 = ax.hist(data,color ='r',rwidth=0.2)
        fig.savefig(outputPath+'.png')   
class simulateFCM:
    def createFCM(self,dimensions,residue_Closeto_activeSite):
        """This function creates a random unbound and bound FCM. The residues close to 
            active site have the maximum number of 1's """
        b = np.random.rand(dimensions,dimensions)
        boundFCM = (b + b.T)/2
        unBoundFCM = np.empty_like (boundFCM)
        unBoundFCM[:] = boundFCM
        maximum = 0
        max_list = []
        for x in range(len(boundFCM[0])):
            max_list.append(sum(i>0.85 for i in boundFCM[x]))
        maximum = max(max_list)
        print maximum
        for x in residue_Closeto_activeSite:
            for x2 in range(maximum):
                boundFCM[x][x2] = 1
        
        return boundFCM,unBoundFCM
        
    def createFCM_2(self,dimensions,residue_Closeto_activeSite):
        """ This is another way of simulating """
        b = np.random.rand(dimensions,dimensions)
        b = np.zeros
        boundFCM = (b + b.T)/2
        unBoundFCM = np.empty_like (boundFCM)
        unBoundFCM[:] = boundFCM
        permuatations = itertools.permutations(residue_Closeto_activeSite,2)
        for x in permuatations:
            boundFCM[x[0]][x[1]] = 0.5
        return boundFCM,unBoundFCM
        
class MyPlots:
    def simplescatterPlot(self, xValues, yValues, xlabel,ylabel,title,outputPath,count):
        #xticks = np.arange(min(xValues),max(xValues)+1,1)
        #yticks = np.arange(min(yValues),max(yValues)+1,0.5)
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle(title,fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel, fontsize= 20)
        ax.set_ylabel(ylabel, fontsize =20)
        #ax.set_xticks(xticks)
        #ax.set_yticks(yticks)
        ax.set_xlim(min(xValues),max(xValues)+1)
        ax.set_ylim(min(yValues),max(yValues)+1)
        ax.tick_params(axis = 'x', labelsize=30)
        ax.tick_params(axis = 'y', labelsize=30)
        pt = ax.scatter(xValues,yValues,color='r')      
        fig.savefig(outputPath+'.png')          
        
        
    def simpleBarPlot(self, xaxis, yValues, xlabel,ylabel,xticks,title,xticklables,outputPath,count):
        xticks = np.arange(1,6,1)
        fig = plt.figure(count, figsize=(30, 20))
        fig.suptitle(title,fontsize= 30)
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel, fontsize= 20)
        ax.set_ylabel(ylabel, fontsize =20)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklables,rotation = 'vertical')
        ax.tick_params(axis = 'x', labelsize=20)
        ax.tick_params(axis = 'y', labelsize=30)
        pt = ax.bar(xaxis,yValues,width=0.05,color='r')      
        fig.savefig(outputPath+'.png')     
        
class statisticalTests():
    def Specificity(self,list1, list2):
        '''
        Compute specificity from two edgelists (list of edges between nodes/labes)
    
        :param ref: Reference edgelist 
        :type ref: list
        :param test: Test edgelist
        :type test: list
        '''
        # Size of intersection
        A = len(set(list1).intersection(list2))
    
        return float(A) / len(list2)
    
    def Sensitivity(self,list1, list2):
        '''
        Compute sensitivity from two edgelists (list of edges between nodes/labes)
    
        :param ref: Reference edgelist
        :type ref: list
        :param test: Test edgelist
        :type test: list
        '''
        # Size of intersection
        A = len(set(list1).intersection(list2))
    
        return float(A) / len(list1)
    
    def Fscore(self,list1, list2):
        ''' 
        Compute the F-score 
    
        :param sp: Specificity value
        :type sp: float
        :param sn: Sensitivity value
        :type sn: float
        '''
        sp = self.Specificity(list1, list2)
        sn = self.Sensitivity(list1,list2)
        fscore = 2*sp*sn/(sp+sn) 
        return sp,sn, fscore   
