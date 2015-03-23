'''
Created on Mar 11, 2015

@author: Andrius Radvilas Bubelis
'''
from types import DictionaryType

if __name__ == '__main__':
    pass

print 'Hi Andrius!'


from Bio import AlignIO, SeqIO, Phylo
from Bio.Align  import MultipleSeqAlignment
from Clade import Clade
import re
import pdb




'''
input vars
'''
consensus = 1  # 0-1% cutoff per column
min_conserved_length = 10 
treefile  = 'dengue.tre' #'C:\Users\Radvilas\SkyDrive\Documents\SemiTreeData\dengue.tre'
fastafile = 'dengue.fa' #'C:\Users\Radvilas\SkyDrive\Documents\SemiTreeData\dengue.tre'
odfile    = ''
phfile    = ''

#import logging
#logging.basicConfig(level=logging.DEBUG)

cladelist  = list();
cladedict  = dict();
cladenames = list();
fasta      = dict() #temp dict
allfasta   = dict() #fasta dict for all selected clades

#find clades and copy into  myclade
tree = Phylo.read(open(treefile), 'newick')
#read tree's fasta file to seq
seq = AlignIO.read( open(fastafile),"fasta" )


seq_dic = SeqIO.to_dict(seq)

MSAeverycalade = MultipleSeqAlignment([])
nonterminals = tree.get_nonterminals()
for clnode in nonterminals:
    cladename  = str(clnode)
    #s.startswith("CLADE"")
    if re.match("CLADE_*", cladename):
        
        
        fasta.clear() #temp fasta list
           
        
        s = list()
        del s[:] 
        MSA = MultipleSeqAlignment(s) #temp MSA filr

            
        
        #get fasta for each clade
        for cname in clnode.get_terminals():
            name = str(cname) 
            fa   = str(seq_dic[name].seq) # one fasta sequence for given name e.g. DNG3 
            MSA.append(seq_dic[name])  # old method: MSA.add_sequence(name, fa)
            fasta[name]    = fa
            allfasta[name] = fa
            
        #print "xxx" , name, fa[1000:]
            
        #print  MSA[:, len(MSA[0])-20:]
            
        myclade = Clade(cladename,clnode,fasta,MSA)

        myclade.calculate_conservation()

        myclade.find_conserved_regions(min_conserved_length,consensus)

        
        #print 'xxxENDxxx'
        #pdb.set_trace()

        
        ################################
        cladenames.append(cladename)
        cladelist.append(myclade)
        cladedict[cladename] = myclade  
        ################################
        
        myclade.write_clade_to_files() # defaulf: cladename.tre, clade_name.fa      

        MSAeverycalade.extend(MSA)
        
        

all_clades = Clade("ALL_CLADES",tree,allfasta,MSAeverycalade) #to calculate global conservation

all_clades.calculate_conservation()

all_clades.find_conserved_regions(min_conserved_length,consensus)

all_clades.write_clade_to_files()        

print 'Done!'

#print 'xxxENDxxx'
#pdb.set_trace()



