'''
Created on Mar 21, 2015

@author: Andrius Radvilas Bubelis
'''
from Bio import Phylo, AlignIO, SeqIO
from Bio.Align  import MultipleSeqAlignment
import copy
import pdb
from telnetlib import theNULL
from Bio.Phylo import Consensus

class ConsrvBlock(object):

    nID = 0         # As Integer '+= 1  debug , count every conserved position recorded (number of total conserved columns)
    Length = 0      # As Integer 'how many conserved AA in block?
    strAA = ''      # As String 'conserved region AA seq (string)
    startPos = -1     # 'start position of the Block


    def __init__(self, nID=0, startPos=-1, strAA=""):
        '''
        Constructor
        '''
        self.nID =nID
        self.startPos = startPos
        self.strAA = strAA
        self.Length = len(strAA)



class Clade(object):
    '''
    classdocs
    '''
    name = ''
    fasta = dict();
    ConsrvRegion = list(ConsrvBlock());
    
    MSAnogaps = MultipleSeqAlignment([]) #stores no-caped version of MSA
    
    '''Aggregates'''
    frequencies = list() 
    mostfreqaa  = list()
    mostfreqsum = list()

    
    
    
    
    
    def __init__(self, name, tree, fasta, MSA):
        '''
        Constructor
        '''
        self.name      = name
        self.tree      = tree
        self.fasta     = fasta
        self.MSA       = MSA
        self.remove_gaps() # creates self.MSAnogaps object

    def remove_gaps(self):
        #remove columns with all gaps
        self.MSAnogaps = copy.copy(self.MSA)
                
        EOL = len(self.MSAnogaps[0])-1
        for col in range(EOL,0,-1) :
             
            column = self.MSAnogaps[:, col]
            #print  col, column
            
            if self.is_gap_colum(column):
                #remove gap column
                if (col==len(self.MSAnogaps[0])-1):
                    print 1, self.MSAnogaps[:, :col-1]
                    self.MSAnogaps = self.MSAnogaps[:, :col-1]
                    pdb.set_trace()
                    
                elif (col==0):
                    print 2, self.MSAnogaps[:, :]
                    self.MSAnogaps = self.MSAnogaps[:, 1:]
                    pdb.set_trace()
                
                else:
                    self.MSAnogaps = self.MSAnogaps[:, :col-1] + self.MSAnogaps[:, col+1:]

    
    def draw_tree(self):
        print self.name
        Phylo.draw_ascii(self.tree)
        
    def print_fasta(self, seq_name = ''):
        if seq_name == '':
            for sec in self.fasta:
                print self.fasta[sec]
        else:
            print seq_name, self.fasta[seq_name]
    
    def is_gap_colum(self, column):
        for s in column:
            #print s
            if s <> "-":
                return False
        return True
      
        
    def write_clade_to_files(self):
        
        treFile   = self.name + '.tre'
        fasFile   = self.name + '.fa'
        fasFileNG = self.name + '_nogaps' + '.fa'
        
        print self.name, "writing tree to file", treFile
        Phylo.write(self.tree, treFile, 'newick')
        
        print self.name, "writing fasta to file", fasFile
        output_handle = open(fasFile, "w")
        AlignIO.write(self.MSA, output_handle, "fasta")
               
        print self.name, "writing fasta with no gaps to file", fasFileNG
        output_handle2 = open(fasFileNG, "w")
        AlignIO.write(self.MSAnogaps, output_handle2, "fasta")
    
    def calculate_conservation(self):

        self.frequencies = []
        self.mostfreqaa  = []
        
        matrix = [list(self.fasta[key]) for key in self.fasta] # for each fasta seq creates list of amino acids (the list of lists) 
        
        # for each column
        for index,item in enumerate(matrix[0]): 
             
            column = [row[index] for row in matrix] # make as list from each row's index e.g. column[1] will store List of all amino acids at position 1
            
            most_amino = '-'
            maxaacount = 0
            aminocount = {}
            # for each row in one column
            for cell in column:
                cell = str(cell)
                                       
                if cell == "-":
                    aminocount[cell] = 0
                elif cell not in aminocount:
                    aminocount[cell] = 1
                else:
                    aminocount[cell] += 1

                #if new amino acid is more frequent then remember it 
                if maxaacount < aminocount[cell] and cell <> "-": 
                    maxaacount = aminocount[cell]
                    most_amino = str(cell)
                
            maxcount = max([aminocount[key] for key in aminocount])
                           
            #if maxcount == 1:
            #    conservation_score = 0
            #else:
            
            conservation_score = float(maxcount) / len(column)
            
            self.frequencies.append(conservation_score)
            self.mostfreqaa.append (most_amino)
            print index, column, conservation_score, most_amino      


    def find_conserved_regions(self, min_len, consensus):
        
        for f in frequencies()
            
        pass
        #print self.mostfreqaa

'''        
class ConsrvBlock(object):

    nID = 0         # As Integer '+= 1  debug , count every conserved position recorded (number of total conserved columns)
    Length = 0      # As Integer 'how many conserved AA in block?
    strAA = ''      # As String 'conserved region AA seq (string)
    startPos = -1     # 'start position of the Block


    def __init__(self, nID=0, startPos=-1, strAA=""):


'''         
            
            
            
            
    
            
            

        
        
        

         