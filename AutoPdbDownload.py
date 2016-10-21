#This program takes a text file as input whic is downloaded from NCBI blast search. The text file contains the pdb codes
# for all the homologs. This programs automatically downloads all the pdb files of the homologs. It parses the tect file and then downloads 
#uisng the variable url.

#Author: Qazi Zaahirah
#email: qazi.zaahirah50@gmail.com
from sys import argv 
import urllib,os

#name of the text file which has all the homologs
#The same file will hold all the pdb files
fileName = '/media/qazi/MyDrive/Downloads/Homologs.txt'
with open(fileName) as textFile:
        lines = textFile.readlines()
        for i, line in enumerate(lines):
                content = line.split()
                #print content
                if content:
                        myWord = content[0]
                        #print myWord
                        if (myWord.startswith('pdb')):
                                checkWord = myWord.split('|')
                                var = checkWord[1]
                                os.system('wget http://files.rcsb.org/download/%s.pdb'%var)
                                
                        
                
               
                
                 
                 
                                         
                 
        


 
 
 
 
