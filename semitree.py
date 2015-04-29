'''
Created on Mar 11, 2015

@author: Andrius Radvilas Bubelis
'''
from types import DictionaryType
from _ssl import txt2obj
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align  import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Clade import *
from Clade import Clade

import re
import argparse
import getpass
import csv # import order disorder file

if __name__ == '__main__':
    pass

print 'Hello ' + getpass.getuser() + '!'

'''
INPUT ARGS

consensus = 0.9  # 0.0-1.0 (%) cutoff per column
min_conserved_length = 4

treefile            = 'p53all.tre' 
fastafile           = 'p53all.fa'  
order_disorder_file = 'p53all.odo'

treefile            = 'F1QF54.tre' 
fastafile           = 'F1QF54.fa'  
order_disorder_file = "F1QF54.odo"
'''


parser=argparse.ArgumentParser(description="This is a SemiTee pipeline")

parser.add_argument("-t", "--input_tree",       default='F1QF54.tre',    help="newick tree file")
parser.add_argument("-a", "--input_alignment",  default='F1QF54.fa',     help="fasta msa file for the tree")
parser.add_argument("-o", "--input_odo",        default='F1QF54.odo',    help="order/disorder file (or any other discreet property file")
parser.add_argument("-c","--consensus",         default=0.9, type=float, help="consensus % for conservation (0.0-1.0)")
parser.add_argument("-m","--minimum_aa",        default=4,   type=int,   help="the minimum number of amino acids for conserved region")
parser.add_argument("-e","--entropy_gap_weight",default=0.5, type=float, help="entropy gap weight range[0.0-1.0] e.g 0.5 => treat gaps like 50% conserved")


args = parser.parse_args()

treefile             = args.input_tree
fastafile            = args.input_alignment
order_disorder_file  = args.input_odo
consensus            = args.consensus
min_conserved_length = args.minimum_aa
entropy_gap_weight   = args.entropy_gap_weight

print "--------------------------------------"
print "----- INPUT ARGUMENTS ----------------"
print "--------------------------------------"
print " treefile             = " , treefile 
print " fastafile            = ", fastafile 
print " order_disorder_file  = ", order_disorder_file 
print " consensus            = ", consensus 
print " min_conserved_length = ", min_conserved_length 
print " entropy_gap_weight   = ", entropy_gap_weight
print "--------------------------------------"


fileheader = treefile.split('.')[0]

   
odfile    = ''
phfile    = ''


#cladedict    = dict()
cladelist    = list()
fasta        = dict() 
allfasta     = dict() #fasta dict for all selected clades
allodo_dic   = dict() #order disorder dictionary
 


#find clades and copy into  myclade
tree = Phylo.read(open(treefile), 'newick')
#read tree's fasta file to seq
seq = AlignIO.read( open(fastafile),"fasta" )
seq_dic = SeqIO.to_dict(seq)

    
if order_disorder_file <> "":
    # round order/disorder from float to 0 or 1
 
    line_nr = 0
    with open(order_disorder_file) as tsv:
        for line in csv.reader(tsv,  delimiter="\t"):
            line_nr += 1 
            seqname = line[0] 
            if seqname <> "LABELS": #ignore header
                bit_line = ''
                
                pos = 0 
                for c in line:
                    if pos == 0:
                        pos +=1
                    else:
                        #print type(c), c, is_number(c)
                        if is_number(c):
                            bit_line +=  str( int(round(float(c))) )
                        else:
                            bit_line += '-'
                        
                
            
                allodo_dic[seqname] = bit_line # adding to all odo to dict, so later we can crop 
            
    tsv.close()


MSA_everycalade    = MultipleSeqAlignment([])
MSAodo_everycalade = MultipleSeqAlignment([])

nonterminals = tree.get_nonterminals()
for clnode in nonterminals:
    cladename  = str(clnode.name)
    #s.startswith("CLADE"")
    
    if re.match("GROUP_*", cladename, re.IGNORECASE) or re.match("CLADE_*", cladename, re.IGNORECASE) :
        #print "CLADE NAME: " + cladename
            
        fasta.clear() #temp fasta list

        tmp1 = []
        tmp2 = []
        del tmp1[:]
        del tmp2[:] 
         
        MSA = MultipleSeqAlignment(tmp1) #temp MSA
        ODO = MultipleSeqAlignment(tmp2) #temp MSA
        
        odo_list = list()
        del odo_list[:]
        
        #get fasta for each clade
        n = 0
        for cname in clnode.get_terminals():
            
            name = str(cname.name)   
            fa   = str(seq_dic[name].seq)
             
            if len(allodo_dic) > 0:  
                odo_record = SeqRecord(Seq(allodo_dic[name]), name=name)
                n +=1
                ODO.append(odo_record)
              
            MSA.append(seq_dic[name])  # old method: MSA.add_sequence(name, fa)

            fasta[name]    = fa
            allfasta[name] = fa

        ### CREATING CLADE ############################################
        myclade = Clade(fileheader, cladename, clnode, fasta, MSA, ODO, entropy_gap_weight)
        ###############################################################
 
               
        '''Oder/DisOrder'''
        myclade.calculate_odo_conservation() #find frequency for each column 
            
        ''' AA conservation '''
        myclade.calculate_aa_conservation() #find frequency for each column

        myclade.find_aa_conserved_regions(min_conserved_length,consensus)
        
        
        ''' AA entropy '''
        myclade.calculate_aa_entropy()
        
        ''' PhyChem '''
        myclade.calculate_phchem_conservation() #find frequency for each column

        myclade.find_phychm_conserved_regions(min_conserved_length,consensus)
        
        ''' Rate4Site '''
        myclade.run_rate4site() 
            

        ''' OutPut '''
        myclade.print_clade_analysis() # display  

        myclade.write_clade_to_files() # default: cladename.tre, clade_name.fa      

        
        
        #cladedict[cladename] = myclade
        cladelist.append(myclade)
        
        
        MSA_everycalade.extend(MSA)
        MSAodo_everycalade.extend(ODO)

''' +++ ALL_CLADES +++ 
ALL SELECTED CLADES IN ONE (it might be the case when not all clades are seleced in tree)
'''

all_clades = Clade(fileheader, "ALL_CLADES", tree, allfasta, MSA_everycalade, MSAodo_everycalade, entropy_gap_weight)


''' Oder/DisOrder '''
all_clades.calculate_odo_conservation() #find O/D frequency for each column 
    
''' AA '''
all_clades.calculate_aa_conservation() #findAA  frequency for each column

all_clades.find_aa_conserved_regions(min_conserved_length,consensus)

''' AA entropy '''
all_clades.calculate_aa_entropy()

''' PhyChem '''

all_clades.calculate_phchem_conservation() #fin PhyChemd frequency for each column

all_clades.find_phychm_conserved_regions(min_conserved_length,consensus)

''' Rate4Site '''
all_clades.run_rate4site()

''' OutPut '''
all_clades.print_clade_analysis() # display  

all_clades.write_clade_to_files() #       


#write rate4site for all clades in one file
txt = ""
for clade in cladelist:
    txt +=  clade.name + ' alpha = ' + str(clade.Alpha) + '\n' 
 
    
txt += 'Pos \t'
clade_names = ''
for clade in cladelist:
    clade_names += clade.name + '\t'
      
txt += clade_names
txt += '\n' 
    
for c in range(0,len(MSA_everycalade[0])):
    txt += str(c+1)
    line = ""
    
    for clade in cladelist:
        #print clade.name, len(clade.rate4site_list)
        
        sSCORE = str(clade.rate4site_list[c].SCORE) 
        txt += '\t' + sSCORE
        line += '\t' + sSCORE 
    #print c, line 
    txt += '\n'   

print "\nFILES CREATED:"
print "Rate4Site global summary: " + fileheader + "_r4s.txt"

f = open(fileheader + "_r4s.txt", "w")
f.write(txt)
f.close()


for clade in cladelist:
    print clade.name
    clade.print_files_created()


print 'Done!'


#print 'xxxENDxxx'
#pdb.set_trace()


