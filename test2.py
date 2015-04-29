'''
Created on Mar 29, 2015
test 
@author: Radvilas

from numpy import *

a = array([1,2,7,5,6])
alfa = a.mean()**2/a.var()

print alfa


list1 = ['1', '-', '0']
str1 = ''.join(list1)
print str1


import collections
a = ['-','-','-','-','-'] #,'T', 'Y','Y']

counter=collections.Counter(a)
print(counter)
# Counter({1: 4, 2: 4, 3: 2, 5: 2, 4: 1})
print(counter.values())
# [4, 4, 2, 1, 2]
print counter.items(), "keys"
# [1, 2, 3, 4, 5]

print(counter.most_common(), "+++", len(counter.keys()))
# [(1, 4), (2, 4), (3, 2)]

dic = dict(counter.items())

print dic['-']

for aa in counter :
    print aa, counter[aa]
    
s = list(counter.most_common())[0][0]

print s
'''
import math

len_column = 10
len_column_nogaps = 1  #AA
len_column_gaps = len_column - len_column_nogaps
  
FrGap =  float(len_column_gaps) / len_column
FrAA  =  1 - FrGap 


d = int(math.ceil(float(len_column_gaps) / (1.0/0.5)))

P =  (float(len_column_gaps) / len_column)  /d 
 
print FrGap , P , d, P* d 
e = 0
for n in range(1,d+1):
    e += P * math.log(P, 2)
    
print e*-1


    
    