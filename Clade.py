'''
Created on Mar 21, 2015

@author: Andrius Radvilas Bubelis
'''
from Bio import Phylo, AlignIO, SeqIO
from Bio.Align  import MultipleSeqAlignment
import copy
import pdb
import collections
import numpy as np
from subprocess import os
from collections import Counter
import math

from sys import platform as _platform

def is_number(s):
    try:
        float(s) # for int, long and float
    except ValueError:
        try:
            complex(s) # for complex
        except ValueError:
            return False
    return True



class ConsrvBlock(object):

    nID = 0          # As Integer '+= 1  debug , count every conserved position recorded (number of total conserved columns)
    Length = 0       # As Integer 'how many conserved AA in block?
    strAA = ''       # As String 'conserved region AA seq (string)
    startPos = -1    # 'start position of the Block
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
        self.entropy_gap_weight = 0.5

class R4site(object):
    POS   = 0
    AA   = ''
    SCORE = 0
    STD   = 0
    CNAME = 'n/a'
    
    def __init__(self, CNAME, POS, AA, SCORE, STD):
        self.POS   = POS 
        self.AA    = AA
        self.SCORE = SCORE
        self.STD   = STD
        self.CNAME = CNAME


class Clade(object):
    '''
    classdocs
    '''
    bRate4SiteDidRun = False
    fileheader = ''
    Alpha = 0 #rate4site
    seqlen = 0
    name = ''
    maxseqlen = 0
    longest_seq_name = ''
    minseqlen = -1
    fasta = {} #disc
    fasta_nogaps = {}
    
    #odo_list            = []
    AAConsrvRegions     = [] # list of objects AA conservation
    PhChemConsrvRegions = [] # list of objects PhyChem conservation
    rate4site_list      = [] # list of objects Rate4Site - stores all columns
    
    rowcount = 0
    colcount = 0
    
    MSAnogaps = MultipleSeqAlignment([]) #stores no-caped version of MSA
    
    
    '''Aggregates'''
    
    #Amino acids
    frequencies        = list() #AA
    frequencies_nogaps = list() #AA
    mostfreqaa         = ''     #AA
    aa_entropy            = list()
    aa_negentropy_normalized = list() # normalized by max and reversed 1-e  e.g. 0 - highest , 1 - lowest  , so would be easy visually compare with conservation
    aa_entropy_nogaps     = list()
    #mostfreqsum        = list() #AA

    #PhyChem
    phchfrequencies_list        = list() 
    phchfrequencies_nogaps_list = list()
    phchmostfreq_str            = ''
    #phchmostfreqsum             = list()

    #Order Disorder
    odofrequencies_list        = list()
    odofrequencies_nogaps_list = list() 
    odomostfreq_str            = ''
    odohmostfreqsum_list       = list()
    
    ''''''''''''''''''

    min_len   = 5
    consensus = 1.0  
    
    
    def __init__(self, fileheader, name, tree, fasta, MSA, ODO, entropy_gap_weight):
        '''
        Constructor
        '''
        self.fileheader  = fileheader
         
        self.name   = name
        self.tree   = tree
        self.fasta  = fasta
        self.MSA    = MSA
        self.ODO    = ODO
        
        self.entropy_gap_weight = entropy_gap_weight 
        
        
        #self.do_list              = []
        self.AAConsrvRegions      = []
        self.PhChemConsrvRegions  = []
        self.rate4site_list       = []
        
        
        self.frequencies        = []
        self.frequencies_nogaps = []
        
        self.aa_entropy            = []
        self.aa_entropy_nogaps     = []
        self.aa_negentropy_normalized = []

    
        #PhyChem
        self.phchfrequencies_list        = [] 
        self.phchfrequencies_nogaps_list = []
        #self.phchmostfreqsum             = []
    
        #Order Disorder
        self.odofrequencies_list        = []
        self.odofrequencies_nogaps_list = [] 
        #self.odohmostfreqsum_list       = []        

        self.rowcount  = len(fasta)
        self.colcount  = len(self.MSA[0])
        
        self.remove_gaps() # creates self.MSAnogaps object, and fasta_nogaps dict
    
        self.treFile         = self.fileheader + '_' + self.name + '.tre'
        self.fasFile         = self.fileheader + '_' + self.name + '.fa'
        self.fasFileNG       = self.fileheader + '_' + self.name + '_nogaps' + '.fa'
        self.consrvFile      = self.fileheader + '_' + self.name + '_conserved_reg' + '.txt'
        self.consrvPyChFile  = self.fileheader + '_' + self.name + '_conserved_phychm_reg' + '.txt'
        self.summaryFile     = self.fileheader + '_' + self.name + '_summary.txt'
        self.csvFile         = self.fileheader + '_' + self.name + '.csv'
        self.r4sFile         = self.fileheader + '_' + self.name + ".r4s"

    
    def make_summary(self):
    
        summary =  '\n******** ' + self.fileheader  + ' ** ' + self.name +' CLADE ANALYSIS ***************'
        summary += '\nMIN REGION LEN: ' + str(self.min_len) + '\nCONSENSUS: ' + str(self.consensus) 
        summary += '\nShortest AA Length=' + str(self.minseqlen)
        summary += '\nLongest AA Length=' + str(self.maxseqlen)
        summary += '\nCLADE Total Length=' + str( self.seqlen)
        
        if self.bRate4SiteDidRun:
            summary += '\nRate4Site Alpha=' + str( self.Alpha )
        
        
        ar   = np.array(self.frequencies_nogaps)
        Alpha = ar.mean()**2/ar.var()
        summary += '---------------------------------------------------'
        summary += '\nAmino Acid Conservation, number of conserved regions:' + str( len(self.AAConsrvRegions))
        summary += '\nAmino Acid Conservation, Mean:' + str( ar.mean() )
        summary += '\nAmino Acid Conservation, STD:' + str( ar.std() )
        summary += '\nAmino Acid Conservation, Alpha:' + str( Alpha )
        summary += '\n---------------------------------------------------'

        ar   = np.array(self.aa_entropy_nogaps)
        Alpha = ar.mean()**2/ar.var()
        summary += '\nAmino Acid Entropy, Mean:' + str( ar.mean() )
        summary += '\nAmino Acid Entropy, STD:' + str( ar.std() )
        summary += '\nAmino Acid Entropy, Alpha:' + str( Alpha )
        summary += '\n---------------------------------------------------'
                     
        ar   = np.array(self.phchfrequencies_nogaps_list)
        Alpha = ar.mean()**2/ar.var()
        summary += '\nPhysicochemical Conservation, number of conserved regions:' + str( len(self.PhChemConsrvRegions))
        summary += '\nPhysicochemical Conservation, Mean:' + str( ar.mean() )
        summary += '\nPhysicochemical Conservation, STD:' + str( ar.std() )
        summary += '\nPhysicochemical Conservation, Alpha:' + str( Alpha )
        summary += '\n---------------------------------------------------'
        
        ar   = np.array(self.odofrequencies_nogaps_list)
        Alpha = ar.mean()**2/ar.var()
        summary += '\nO/D Conservation, Mean:' + str( ar.mean() )
        summary += '\nO/D Conservation, STD:' + str( ar.std() )
        summary += '\nO/D Conservation, Alpha:' + str( Alpha )


        self.clade_summary = summary


    def print_clade_analysis(self):
        
        self.make_summary()
        
        print self.clade_summary 

        
    
    def remove_gaps(self):
        

        
        #remove columns with all gaps
        self.MSAnogaps = copy.copy(self.MSA)
        
                
        EOL = len(self.MSAnogaps[0]) - 1
        #print "EOL", EOL
        
        for col in range(EOL,0,-1) :
            
            #print  "col=", col, 
            lastcol = len(self.MSAnogaps[0])-1
            if (lastcol < col):
                col = lastcol
            
            column = self.MSAnogaps[:, col]
            
                 
            if self.is_gap_colum(column):
                #remove gap column
                #print  "is_gap_colum = TRUE",col 
                
                if (col==len(self.MSAnogaps[0])-1):
                    #print "*1*", self.MSAnogaps[:, :col-1]
                    self.MSAnogaps = self.MSAnogaps[:, :col-1]

                    
                elif (col==0):
                    #print "*2*", self.MSAnogaps[:, 1:]
                    self.MSAnogaps = self.MSAnogaps[:, 1:]
                    #print '2xxxENDxxx'
                    #pdb.set_trace()   
                
                else:
                    #print "*3*"
                    self.MSAnogaps = self.MSAnogaps[:, :col-1] + self.MSAnogaps[:, col+1:]
                    
                    
                    
        #remove columns with all gaps for ODO=============================================
        self.ODOnogaps = copy.copy(self.ODO)
        
                
        EOL = len(self.ODOnogaps[0]) - 1
        #print "EOL", EOL
        
        for col in range(EOL,0,-1) :
            
            #print  "col=", col, 
            lastcol = len(self.ODOnogaps[0])-1
            if (lastcol < col):
                col = lastcol
            
            column = self.ODOnogaps[:, col]
            
                 
            if self.is_gap_colum(column):
                #remove gap column
                #print  "is_gap_colum = TRUE",col 
                
                if (col==len(self.ODOnogaps[0])-1):
                    #print "*1*", self.MSAnogaps[:, :col-1]
                    self.ODOnogaps = self.ODOnogaps[:, :col-1]

                    
                elif (col==0):
                    #print "*2*", self.MSAnogaps[:, 1:]
                    self.ODOnogaps = self.ODOnogaps[:, 1:]
                    #print '2xxxENDxxx'
                    #pdb.set_trace()   
                
                else:
                    #print "*3*"
                    self.ODOnogaps = self.ODOnogaps[:, :col-1] + self.ODOnogaps[:, col+1:]               
                    
        
        # make dict fasta_nogaps==============================================
        self.maxseqlen = 0
        self.minseqlen = -1
        self.seqlen = len(self.MSAnogaps[0])
        for seq in self.MSAnogaps:
            s = str(seq.seq)
            self.fasta_nogaps[seq.name] = s
            seqnogaps = s.replace('-','')
            l = len(seqnogaps)

            if self.maxseqlen < l:
                self.maxseqlen = l
                self.longest_seq_name = seq.name 
                  
            if self.minseqlen < 0 or l < self.minseqlen: 
                self.minseqlen = l
            

    #remove_gaps() end#########################
    
        
    def draw_tree(self):
        #print self.name
        Phylo.draw_ascii(self.tree)
        
    def print_fasta(self, seq_name = ''):
        if seq_name == '':
            for sec in self.fasta:
                print self.fasta[sec]
        else:
            print seq_name, self.fasta[seq_name]
    
    def is_gap_colum(self, column):
        for s in column:
            if s <> "-":
                return False
        return True
      
    def parse_rate4site_to_list(self,r4sfile):
        POS   = 0
        n = 0    

        if os.path.isfile(r4sfile) == False: 
            print "!!!!!Rate4site output file not found: " + self.r4sFile
            
        with open(r4sfile) as f:

            for fline in f:
                
                n +=1
                
                line = str(fline)
                
                if len(line) > 0:
                    if '#The alpha parameter' in line:
                        self.Alpha = line[20:]
                        #print "r4s.Alpha = ", self.Alpha , n
                    
                    table = line.split()
                    #print len(table) , "+++++++++"
                    if len(table) > 0:
                        if is_number(table[0]):
                            
                            POS   = int(table[0])
                            AA    = table[1]
                            SCORE = float(table[2])
                            try:
                                STD   = float(table[5])
                            except:
                                try:
                                    STD   = float(table[4])
                                except:
                                    STD   = float(table[6])
                            
                            
                            s =  str(POS) +  '+' +  AA + '+' + str(SCORE) +  '+' + str(STD)
                            
                            self.rate4site_list.append( R4site(self.name, POS, AA, SCORE, STD) )
                            #print self.name, s, len(self.rate4site_list)
                            
                
        f.close

        if POS < len(self.MSAnogaps[0]):
            for p in range(POS+1, len(self.MSAnogaps[0])+1):
                self.rate4site_list.append( R4site(p, 'X', 0.777, 0.25) )
                          
    def run_rate4site(self, bRunRate4Site = True):
        
        if bRunRate4Site:
            
            if "linux" in _platform:
                cmdtxt = "rate4site -t " + self.treFile + " -s " + self.fasFile + " -o " + self.r4sFile
                print cmdtxt
                os.system(cmdtxt)
                #cmdtxt = "unix2dos.exe " + r4sFile
                #print cmdtxt
            
            if os.path.isfile(self.r4sFile): 
                
                self.parse_rate4site_to_list (self.r4sFile)
            
                self.bRate4SiteDidRun = True
            else:
                print "######Rate4site output file not found: " + self.r4sFile
                print "######Rate4site output file not found: " + self.r4sFile
                print "######Rate4site output file not found: " + self.r4sFile
                
                
            
    def print_files_created(self):
        
        print self.name, "tree: ", self.treFile
        
        print self.name, "fasta: ", self.fasFile
               
        print self.name, "fasta without gaps: ", self.fasFileNG
        
        print self.name, "amino acid conserved region summary:", self.consrvFile

        print self.name, "physicochemical conserved region summary: ", self.consrvPyChFile

        print self.name, "clade summary: ", self.consrvFile

        print self.name, "clade conservation CSV file: ", self.csvFile


    def write_clade_to_files(self, bRunRate4Site = True):
        
        if self.bRate4SiteDidRun == False:
            bRunRate4Site = False
        
        #print self.name, "writing tree to file", self.treFile
        Phylo.write(self.tree, self.treFile, 'newick')
        
        #print self.name, "writing fasta to file", self.fasFile
        output_handle = open(self.fasFile, "w")
        AlignIO.write(self.MSA, output_handle, "fasta")
               
        #print self.name, "writing fasta with no gaps to file", self.fasFileNG
        output_handle2 = open(self.fasFileNG, "w")
        AlignIO.write(self.MSAnogaps, output_handle2, "fasta")
        
        #print self.name, "writing amino acid conserved region summary", self.consrvFile
        f = open(self.consrvFile, "w")
        line = "CLADE NAME: " + self.name + "\nMIN REGION LEN: " + str(self.min_len) + "\nCONSENSUS: " + str(self.consensus) + "\nREGION COUNT:" + str(len(self.AAConsrvRegions)) + "\n\n"
        f.write(line)
        line = "ID\t Pos\t Len\t AA\n"
        f.write(line)
        for reg in self.AAConsrvRegions:
            line = str(reg.nID) + '\t ' + str(reg.startPos+1)  + '\t ' +  str(reg.Length)  + '\t ' +  reg.strAA + '\n'
            f.write(line)
        f.close()

        #print self.name, "writing physicochemical conserved region summary", self.consrvPyChFile
        f = open(self.consrvPyChFile, "w")
        line = "CLADE NAME: " + self.name + "\nMIN REGION LEN: " + str(self.min_len) + "\nCONSENSUS: " + str(self.consensus) + "\n\n"
        f.write(line)
        line = "ID\t Pos\t Len\t PhysicoChemical property\n"
        f.write(line)
        for reg in self.PhChemConsrvRegions:
            line = str(reg.nID) + '\t ' + str(reg.startPos+1)  + '\t ' +  str(reg.Length)  + '\t ' +  reg.strAA + '\n'
            f.write(line)
        f.close()

        #print self.name, "writing clade summary", self.consrvFile
        f = open(self.summaryFile, "w")
        f.write(self.clade_summary)
        f.close()

        EOL = len(self.frequencies) -1
        hd = 'COLUMNS:,'
        if bRunRate4Site:
            f0 = 'RATE4SITE:,'
        
        f1 = 'AA_FREQ:,'
        f2 = 'PHCEM_FREQ:,'
        f3 = 'ODO_FREQ:,'
        f4 = 'AA_NEGENTROPY_normalized:,'
        f5 = 'AA_ENTROPY_Log2:,'
        
        
        for col in range(0,EOL-1):
            hd += str(col) + ','
            # if column is emplty
            if self.odofrequencies_list[col] == 0:
                if bRunRate4Site:
                    f0 += ','
                f1 += ','
                f2 += ','
                f3 += ','
                f4 += ','
                f5 += ','
            else:
                if bRunRate4Site:
                    f0 += str(self.rate4site_list[col].SCORE) + ','
                f1 += str(self.frequencies[col]) + ','
                f2 += str(self.phchfrequencies_list[col]) + ','
                f3 += str(self.odofrequencies_list[col]) + ','
                f4 += str(self.aa_negentropy_normalized[col]) + ','
                f5 += str(self.aa_entropy[col]) + ','

        hd += str(EOL) + '\n'                
        if self.odofrequencies_list[EOL] == 0:
            if bRunRate4Site:
                f0 += '\n'
           
            f1 += '\n' 
            f2 += '\n'
            f3 += '\n'
            f4 += '\n'
            f5 += '\n'
            
        else:
            if bRunRate4Site:
                f0 += str(self.rate4site_list[EOL].SCORE) + '\n'
            f1 += str(self.frequencies[EOL]) + '\n' 
            f2 += str(self.phchfrequencies_list[EOL]) + '\n'
            f3 += str(self.odofrequencies_list[EOL]) + '\n'
            f4 += str(self.aa_negentropy_normalized[EOL]) + '\n'
            f5 += str(self.aa_entropy[EOL]) + '\n'
        
        if bRunRate4Site:
            txt = hd + f0 + f1 + f2 + f3 + f4 + f5
        else:
            txt = hd + f1 + f2 + f3 + f4 + f5
            

        #print self.name, "writing clade conservation CSV file", self.csvFile
        f = open(self.csvFile, "w")
        f.write(txt)
        f.close()
        
 
    def ppptest(self):          
        if self.name == 'CLADE_WNV':
            for reg in self.AAConsrvRegions:
                line = str(reg.nID) + '\t ' + str(reg.startPos+1)  + '\t ' +  str(reg.Length)  + '\t ' +  reg.strAA + '\n'
                print line
            
      
    def calculate_aa_conservation(self):

        self.frequencies = []
        self.frequencies_nogaps = []
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
                           
            #if column is emplty
            if most_amino == '-':
                self.frequencies.append(1)
            else:
                conservation_score = float(maxcount) / len(column)
                self.frequencies_nogaps.append(conservation_score)
                self.frequencies.append(conservation_score)
                
            
            
            self.mostfreqaa += most_amino
            #if self.name =="ALL_CLADES":
            #    print index, column, conservation_score, most_amino 
            

    def calculate_aa_entropy(self):

        self.aa_entropy        = []
        self.aa_entropy_nogaps = []
        
        matrix = [list(self.fasta[key]) for key in self.fasta] # for each fasta seq creates list of amino acids (the list of lists) 
        max_entropy = 0
        # for each column
        for index,item in enumerate(matrix[0]):
             
            column = [row[index] for row in matrix] # make as list from each row's index e.g. column[1] will store List of all amino acids at position 1

            aa_counts   =  Counter(column)  # frequency of all amino acids
            sum_entropy = 0
              
            len_column_nogaps = 0  
            len_column        = len(column)
            
            for n in column:
                if n <> '-':
                    len_column_nogaps += 1
                    
            if len_column_nogaps == 0: # if empty column 
                
                self.aa_entropy.append('')
                   
                #print self.name, index, '-'
            else:  
                for aa in aa_counts:
                    # let's count entropy for column
                    
                    if (aa <> "-"): # ignore gaps   
                        #    expectation = 0.5  # default expectation for calculating entropy = 0.5 (user can change this weight)
        
                        expectation = float(aa_counts[aa])/len_column #  expectation for each aa = (number of amino acids) / (total column length)
                        
                        sum_entropy += expectation * math.log(expectation,2)
                        
                        #print aa, aa_counts[aa], expectation, expectation * math.log(expectation,2)
                
                len_column_gaps = len_column - len_column_nogaps
  
                if len_column_gaps > 0 :
                    
                    
                    #print len_column_gaps, len_column, (1.0/entropy_gap_weight)
                    
                    gap_partitions = int(math.ceil(float(len_column_gaps) / (1.0/self.entropy_gap_weight)))
    
                    gap_partition_probability =  (float(len_column_gaps) / len_column) /  gap_partitions
    
                    for n in range(1,gap_partitions+1):
                        sum_entropy += gap_partition_probability * math.log(gap_partition_probability, 2)
                     
                if sum_entropy <> 0:
                    sum_entropy = sum_entropy * (-1)
                
                
                if max_entropy < sum_entropy:
                    max_entropy = sum_entropy
                
                self.aa_entropy_nogaps.append(sum_entropy)
                self.aa_entropy.append(str(sum_entropy))
                
        # normalized by max_entropy  
        self.aa_negentropy_normalized = []
        for e in self.aa_entropy:
            if e == '':
                self.aa_negentropy_normalized.append('')
            else:
                if max_entropy == 0:
                    self.aa_negentropy_normalized.append('1')
                else:
                    en = 1 - (float(e) / max_entropy)
                    self.aa_negentropy_normalized.append(str(en))
                     
        
                


              
        #print self.aa_entropy

    def calculate_phchem_conservation(self):

        PhyChmProperties = {
        
            #Hydrophobic
            "A":"F", 
            "M":"F",
            "V":"F",
            "I":"F",
            "L":"F",
            "W":"F",
            "F":"F",
            "Y":"F",

            #Hydroplilic
            "N":"H",  
            "Q":"H",
            "G":"H",
            "C":"H",
            "S":"H",
            "T":"H",
            "P":"H",
             
            #Negativey-charged
            "D":"N", 
            "E":"N",
            
            #Positively-charged
            "H":"P", 
            "R":"P",
            "K":"P",
            
            "-":"-"
            
            }

        
        self.phcmfrequencies_list = []
        self.phchfrequencies_nogaps_list = []
        self.phchmostfreqaa_str  = ''
        
        matrix = [list(self.fasta[key]) for key in self.fasta] # for each fasta seq creates list of amino acids (the list of lists) 

        c=0
        # for each column
        for index,item in enumerate(matrix[0]): 
             
            column = [row[index] for row in matrix] # make as list from each row's index e.g. column[1] will store List of all amino acids at position 1
            
            most_phchem = '-'
            maxphchcount  = 0
            phchcount   = {}
            c +=1 
            # for each row in one column
            for cell in column:
                
                cell   = str(cell)
                if cell in PhyChmProperties:
                    phchem = PhyChmProperties[cell]
                else:
                    phchem = "-"  

                if phchem == "-":
                    phchcount["-"] = 0
                elif phchem not in phchcount:
                    phchcount[phchem] = 1
                else:
                    phchcount[phchem] += 1

                #if new PhCm count is more frequent then remember it 
                if maxphchcount < phchcount[phchem] and phchem <> "-": 
                    maxphchcount = phchcount[phchem] # number of most frequent character
                    most_phchem  = str(phchem) # the most frequent character
                
            maxcount = max( [phchcount[key] for key in phchcount] ) 
                           
            if most_phchem == '-':
                self.phchfrequencies_list.append(1) # for calculating conserved regions
            else:
                conservation_score = float(maxcount) / len(column)
                self.phchfrequencies_list.append(conservation_score)
                self.phchfrequencies_nogaps_list.append(conservation_score)
                 
            
            self.phchmostfreqaa_str += most_phchem

     
    
    def calculate_odo_conservation(self):

        if len(self.ODOnogaps) > 0:
            
            self.odofrequencies_nogaps_list = []
           
            EOL = len(self.ODOnogaps[0]) 
                
            for col in range(0,EOL):
    
                column = self.ODOnogaps[:, col]
                #print self.name, column
                counter=collections.Counter(column)
     
                #print(counter)
    
                counter_dic = dict(counter.items())
    
                try:
                    counter_dic['1'] +=0
                except:
                    counter_dic['1'] =0
                    
                try:
                    counter_dic['0'] +=0
                except:
                    counter_dic['0'] =0
    
                
                if counter_dic['1'] + counter_dic['0'] == 0:
                    conservation_score = 1
                    #print self.name, '-', conservation_score
                    
                elif counter_dic['1'] < counter_dic['0']:
                    conservation_score = float(counter_dic['0']) / len(column)
                    #print self.name, '0', conservation_score
                    
                else: 
                    conservation_score = float(counter_dic['1']) / len(column)
                    #print self.name, '1', conservation_score
                    
                self.odofrequencies_nogaps_list.append(conservation_score)
            
            # do same (calculate odo frequency) with gaps    
            self.calculate_odo_conservation_with_gaps() 
        


    def calculate_odo_conservation_with_gaps(self):

        if len(self.ODO) > 0:
            
            self.odofrequencies_list = []
            self.odomostfreqaa  = ''
                        
            EOL = len(self.ODO[0])

            for col in range(0,EOL):
                #print col, EOL
                column = self.ODO[:, col]
                #print self.name, column
                counter=collections.Counter(column)
     
                #print(counter)
    
                counter_dic = dict(counter.items())
    
                try:
                    counter_dic['1'] +=0
                except:
                    counter_dic['1'] =0
                    
                try:
                    counter_dic['0'] +=0
                except:
                    counter_dic['0'] =0
    
                
                if counter_dic['1'] + counter_dic['0'] == 0:
                    self.odomostfreqaa += '-'
                    conservation_score = 0
                    #print self.name, '-', conservation_score
                    
                elif counter_dic['1'] < counter_dic['0']:
                    self.odomostfreqaa += '0'
                    conservation_score = float(counter_dic['0']) / len(column)
                    #print self.name, '0', conservation_score
                    
                else: 
                    self.odomostfreqaa += '1'
                    conservation_score = float(counter_dic['1']) / len(column)
                    #print self.name, '1', conservation_score
                    
                self.odofrequencies_list.append(conservation_score)

            
            #print self.odofrequencies_nogaps_list
            #print self.name, self.odomostfreqaa
              
            

    def find_aa_conserved_regions(self, min_len, consensus):

        nID = 0
        bNewInProgress = False
        #self.AAConsrvRegions = ()
        self.min_len   = min_len
        self.consensus = float(consensus)  
        
        #print self.name
        for pos,frq in enumerate(self.frequencies):
            #print pos,frq,self.mostfreqaa[pos], consensus
            if (frq >= self.consensus): #0.9 > 0.8
                if bNewInProgress == False:
                    
                    bNewInProgress = True
                    startPos = pos
                
            else:
                if bNewInProgress:
                    bNewInProgress = False
                    strAA = self.mostfreqaa[startPos:pos] # from startPos to pos 
                    nogapsRegion = strAA.replace('-','')
                    if (len(nogapsRegion) >= min_len): #if min length of conserved region is long enough then save it to ConsrvRegion list
                        nID += 1
                        self.AAConsrvRegions.append(ConsrvBlock(nID,startPos,strAA,nogapsRegion))
                        #if self.name =="ALL_CLADES":
                        #    print nID,startPos,strAA,nogaps     


    def find_phychm_conserved_regions(self, min_len, consensus):
            
            nID = 0
            bNewInProgress = False
            #self.AAConsrvRegions = ()
            self.min_len   = min_len
            self.consensus = float(consensus)  
            #print self.name
            
            for pos,frq in enumerate(self.phchfrequencies_list):
                #print pos,frq,self.phchfrequencies[pos], consensus
                if (frq >= self.consensus): #0.9 > 0.8
                    if bNewInProgress == False:
                        
                        bNewInProgress = True
                        startPos = pos
                    
                else:
                    if bNewInProgress:
                        bNewInProgress = False
                        strAA = self.phchmostfreqaa_str[startPos:pos]

                        nogapsRegion = strAA.replace('-','')
                        if (len(nogapsRegion) >= min_len): #if min length of conserved region is long enough then save it to ConsrvRegion list
                            nID += 1
                            self.PhChemConsrvRegions.append(ConsrvBlock(nID,startPos,strAA,nogapsRegion))
                            #if self.name =="ALL_CLADES":
                            #    print nID,startPos,strAA,nogaps     

   
            
            
            
    
            
            

        
        
        

         