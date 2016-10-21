
#This program takes the PDB files which have multiple chains and divides them into a single chain
# This program uses PDBnet to read the PDB files
#Author: Qazi Zaahirah
from sys import argv 
import urllib,os
from labblouin import PDBnet #(https://github.com/LabBlouin/LabBlouinTools/tree/master/labblouin)

#this part of the code gets all the chains for a single pdb and forms different PDBs with single chain

os.chdir('/Folder which has the PDB files')
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
        c.WriteAsPDB('Output Folder for the new PDB chain files'+newName)


        

#this part of program reads a text file which has the corresponding chains and then reads the chain no.
# the program then creates the chain into a PDb file

fileName = 'Text file location which has the chains vlaues'
mylist =[]
count = 0
with open(fileName) as textFile:
    lines = textFile.readlines()
    for i, line in enumerate(lines):
        content = line.split()
        if content:
            myfile = content[0]
            myChain = content[1]
            mylist.append([myfile,myChain])
            count =count+1
print mylist

os.chdir('Input Folder for PDB files')
fileNames = os.listdir('.')
filec =0
print len(mylist)
for f in fileNames:
    #print f
    p = PDBnet.PDBstructure(f)
    for j in range(len(mylist)):
        #print mylist[j][0]
        name = mylist[j][0]+'.pdb'
        if name==f:    
            c = p.GetChain(mylist[j][1])
            newname = mylist[j][0]+mylist[j][1]+'.pdb'
            c.WriteAsPDB(Output Folder for the new PDB chain files'+newname)
            filec=filec+1
    
print 'finish' 
print filec

