"""This code downloads xml files. It Writes PDB codes from the file names in the folder file path into a text file and uses the PDBcodes in wget"""

#author: Qazi Zaahirah
#email: qazi.zaahirah50@gmail.com
from sys import argv 
import urllib,os


path = 'Folder which has the PDB Files'
os.chdir(path)
pdbfiles = os.listdir('.')
filepath = 'Folder where you want to store the PDB files'


for files in pdbfiles:
    os.system('wget http://files.rcsb.org/download/%s.xml.gz'%files[:4])
    f = open(filepath,'a')
    f.write(str(files[:5]))
    f.write('\n')
    f.flush()
f.close()        
