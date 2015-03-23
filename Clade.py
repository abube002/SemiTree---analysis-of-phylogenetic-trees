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
from operator import pos
from Canvas import Line

class ConsrvBlock(object):

    nID = 0         # As Integer '+= 1  debug , count every conserved position recorded (number of total conserved columns)
    Length = 0      # As Integer 'how many conserved AA in block?
    strAA = ''      # As String 'conserved region AA seq (string)
    startPos = -1     # 'start position of the Block
    strAAnogaps = ''


    def __init__(self, nID, startPos, strAA, strAAnogaps):
        '''
        Constructor
        '''
        self.nID =nID
        self.startPos = startPos
        self.strAA = strAA
        self.strAAnogaps = strAAnogaps
        self.Length = len(strAAnogaps)



class Clade(object):
    '''
    classdocs
    '''
    name = ''
    fasta = {};
    ConsrvRegions = []
    rowcount = 0
    colcount = 0
    
    MSAnogaps = MultipleSeqAlignment([]) #stores no-caped version of MSA
    
    '''Aggregates'''
    frequencies = list() 
    mostfreqaa  = ''
    mostfreqsum = list()
    ''''''''''''''''''

    min_len   = 5
    consensus = 1.0  
            
    
    
    
    
    def __init__(self, name, tree, fasta, MSA):
        '''
        Constructor
        '''
        self.name      = name
        self.tree      = tree
        self.fasta     = fasta
        self.MSA       = MSA
        
        self.rowcount  = len(fasta)
        self.colcount  = len(self.MSA[0])
        
        self.remove_gaps() # creates self.MSAnogaps object
        
        self.ConsrvRegions = [];

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
                    #print 1, self.MSAnogaps[:, :col-1]
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
        
        treFile    = self.name + '.tre'
        fasFile    = self.name + '.fa'
        fasFileNG  = self.name + '_nogaps' + '.fa'
        consrvFile = self.name + '_conserved_reg' + '.txt'
        
        print self.name, "writing tree to file", treFile
        Phylo.write(self.tree, treFile, 'newick')
        
        print self.name, "writing fasta to file", fasFile
        output_handle = open(fasFile, "w")
        AlignIO.write(self.MSA, output_handle, "fasta")
               
        print self.name, "writing fasta with no gaps to file", fasFileNG
        output_handle2 = open(fasFileNG, "w")
        AlignIO.write(self.MSAnogaps, output_handle2, "fasta")
    
        
        f = open(consrvFile, "w")
        line = "NAME: " + self.name + "\nMIN REGION LEN: " + str(self.min_len) + "\nCONSENSUS: " + str(self.consensus) + "\n\n"
        f.write(line)
        line = "ID\t Pos\t Len\t AA\n"
        f.write(line)
        for reg in self.ConsrvRegions:
            line = str(reg.nID) + '\t ' + str(reg.startPos+1)  + '\t ' +  str(reg.Length)  + '\t ' +  reg.strAA + '\n'
            f.write(line)
            
        f.close()
 
    def ppptest(self):
        if self.name == 'CLADE_WNV':
            for reg in self.ConsrvRegions:
                line = str(reg.nID) + '\t ' + str(reg.startPos+1)  + '\t ' +  str(reg.Length)  + '\t ' +  reg.strAA + '\n'
                print line
            
    
    
    def calculate_conservation(self):

        self.frequencies = []
        self.mostfreqaa  = ''
        
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
                           
            if most_amino == '-':
                conservation_score = 1
            else:
                conservation_score = float(maxcount) / len(column)
            
            self.frequencies.append(conservation_score)
            self.mostfreqaa += most_amino
            #if self.name =="ALL_CLADES":
            #    print index, column, conservation_score, most_amino      

            

    def find_conserved_regions(self, min_len, consensus):
        
        nID = 0
        bNewInProgress = False
        #self.ConsrvRegions = ()
        self.min_len   = min_len
        self.consensus = float(consensus)  
        
        print self.name
        for pos,frq in enumerate(self.frequencies):
            #print pos,frq,self.mostfreqaa[pos], consensus
            if (frq >= self.consensus): #0.9 > 0.8
                if bNewInProgress == False:
                    
                    bNewInProgress = True
                    startPos = pos
                
            else:
                if bNewInProgress:
                    bNewInProgress = False
                    strAA = self.mostfreqaa[startPos:pos]
                    nogaps = strAA.replace('-','')
                    if (len(nogaps) >= min_len): #if min length of conserved region is long enough then save it to ConsrvRegion list
                        nID += 1
                        self.ConsrvRegions.append(ConsrvBlock(nID,startPos,strAA,nogaps))
                        #if self.name =="ALL_CLADES":
                        #    print nID,startPos,strAA,nogaps     


   
            
            
            
    
            
            

        
        
        

         