from labblouin import PDBnet,FASTAnet #(https://github.com/LabBlouin/LabBlouinTools/tree/master/labblouin)
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

#Author: Qazi Zaahirah

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
        
            
        newdict={aa:[] for aa in aalis}
        
        for aa, entropy in newlist:
            
            newdict[aa].append(entropy)
            
        filepath =compPath+'.txt'
        f = open(filepath,'w')
        f.write(str(newdict))
        return newdict
        
    def BoxPlotAggregateEntropyfromtxt(self,path,dictlist):
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
        
    def calcEntropy (self, matrix):
        """This function just takes a matrix as an input and calculates the entropy of each row which corresponds to a residue """
        
        entropylist = []
        for i in range(matrix.shape[0]):
            
            entropylist.append(-(sum([matrix[i][j]*math.log(matrix[i][j],2) for j in range(matrix.shape[1]) if matrix[i][j] > 0])))
            
        return entropylist
    def BoxPlotEntropyForAA(self,aalis,homolis,entropy,chain,compPath):
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
        #newdict={aa:[] for aa in aalis}
        print newdict
        for aa, entropy in newlist:
            
            newdict[aa].append(entropy)
        labels=newdict.keys()
        print newdict
        data=[]
        
        for res in newdict.keys():
            
            mean=np.mean(newdict[res])
            sd=np.std(newdict[res])
            
            #data1=[mean,sd]
            data1=newdict[res]
            data.append(data1)
        plt.boxplot(data)
        plt.boxplot(data,labels=labels)
        plt.ylim(-1,5)
        #another way to get the labels for the x axis
        #x = np.arange(1,20,1)        
        #plt.xticks(x,labels)
        #plt.show()
        # Create a figure instance
        fig = plt.figure(1, figsize=(9, 6))
        
        # Create an axes instance
        ax = fig.add_subplot(111)
        
        # Create the boxplot
        bp = ax.boxplot(data,labels=labels)
        
        # Save the figure
        fig.savefig(compPath+'.png', bbox_inches='tight')
        filepath =compPath+'.txt'
        f = open(filepath,'w')
        f.write(str(newdict))        
        print(newdict)    
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
        entropylist=[]
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
            entropylist.append(-sum(freqlis))
        return entropylist
    def BoxPlotSeqEntropy(self,path,folderlist):

		#This function reads a txt file that consists of the dictionary" 
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
    def SeqEntropyplot(self,path):
        """ This function reads a pdb file, gets plots entropy vs
            PDB sequence plot. This function also calls distligandresidue from the class active site
            to get the index of the active site"""
        pdb=PDBnet.PDBstructure(path)
        c = activeSite()
        activesiteindex = c.distligandresidue(path, (1.814,-7.495,-14.178))
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
        
        
    def SeqEntropyPlotfor2Sequences(self,path1,path2,fastapath):
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
        for i in index1:
            entropy1.insert(i,-0.01)
        for i,j in enumerate(entropy1):
            if '!' in str(j):
                entropy1[i] = float(str(entropy1[i].replace('!','')))
                asIndex1 = i
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
        plt.xticks(x, alignSeq1,rotation='horizontal')
        plt.tick_params(labelsize=6)
        plt.xlabel('Homologous Sequence')
        plt.ylabel('Entropy')
        plt.plot(sma1,'g-')
        plt.annotate('activesiteB', xy=(asIndex1, sma1[asIndex1]), xytext=(activesiteindex1,sma1[asIndex1]+2),
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
        plt.annotate('activesiteUB', xy=(asIndex2, sma2[asIndex2]), xytext=(activesiteindex2,sma2[asIndex2]+1),
                     arrowprops=dict(facecolor='black', shrink=0.05))
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
    
    
    def distFromLigands(self,pdb,compPath):
        """This function assumes that mg and mn are the active sites of the protein. This function access a new function in PDBnet
               def DistanceToligand(self, ligandcord):
                   return distance.dist(self.GetPosition(),ligandcord) """
        m = pdb.orderofmodels[0]
        model=pdb.GetModel(m)
        residues = model.GetResidues()
        
        homoindex=[r.index for r in residues if r.GetAtoms()[0].tempFactor!=-0.01 ]
        entropy=[r.GetAtoms()[0].tempFactor for r in residues if r.index in homoindex]
        enolaseB=(26.634 ,74.962 ,12.870)
        gdB=(-23.716, -8.568, -8.136)
        mrB =(-20.292,60.052 ,55.442)
        mleB =(1.814,-7.495,-14.178)
        enolaseUB=(12.403,-15.797,14.352)
        gdUB=(115.429,56.667,48.841)
        mrUB=(70.232,-5.565,10.634)
        mleUB=(15.682,-14.294,-91.348)
        distance=[]
        activeSiteindex = activeSite.distligandresidue(self, compPath+'.pdb',mleB)
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
        tup=[[d,e] for d,e in zip(distance,entropy) if d<15]
        randomtup = []
        x = [ a[0] for a in tup]
        y = [ a[1] for a in tup]        
        for i in range(0,len(x)):
            randomtup.append(random.choice(Proteintuple))        
        f = open(compPath+'random.csv', 'w')
        writer = csv.writer(f,delimiter='\t') 
        for t in randomtup:
            writer.writerow(t)
        f.close()
        #np.save(compPath, tup)
        plt.title('Active Site for GlucarateDehydratase Bound')
        plt.ylim(0,8)
        plt.xlabel('DistanceFromActiveSite')
        plt.ylabel('Entropy')

        
        for i in range(0,len(x)):
            randomtup.append(random.choice(Proteintuple))
        rx = [ a[0] for a in randomtup]
        ry = [ a[1] for a in randomtup]
        plt.plot(x,y,"o")
        #plt.show()        
        
    def distligandresidue(self, path,ligandcoordinates): 
        """This function takes a PDB file and ligand coordinates as input and finds the index of the residue 
        which is closest to the ligand"""
        pdb=PDBnet.PDBstructure(path)
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
        print 'hello'
        
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
