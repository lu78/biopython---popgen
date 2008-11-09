#!/usr/bin/env python


"""
>>> import fastPhaseOutputIO
>>> from StringIO import StringIO

# define a short function to convert fastPhase files into fasta:
>>> def convert2fasta(inputfile):
...     for seqrecord in fastPhaseOutputIO.fastphaseoutputIterator(inputfile):
...         print seqrecord.format("fasta")

##### Short fastPhase output file
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T T T T T G A A A C C A A A G A C G C T G C G T C A G C C T G C A A T C T G
... T T T T T G C C C C C A A A A G C G C G T C G T C A G T C T A A G A C C T A
... Ind2  # subpop. label: 6  (internally 1)
... C T T T T G C C C T C A A A A G T G C T G T G C C A G T C T A C G G C C T G
... T T T T T G A A A C C A A A G A C G C T T C G T C A G T A T A C G A T C T A
... END GENOTYPES
... ''')  

>>> convert2fasta(phasefile)
>Ind1_all1
TTTTTGAAACCAAAGACGCTGCGTCAGCCTGCAATCTG
>Ind1_all2
TTTTTGCCCCCAAAAGCGCGTCGTCAGTCTAAGACCTA
>Ind2_all1
CTTTTGCCCTCAAAAGTGCTGTGCCAGTCTACGGCCTG
>Ind2_all2
TTTTTGAAACCAAAGACGCTTCGTCAGTATACGATCTA
<BLANKLINE>

##### fastPhase file with a single marker and 2 individuals
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T
... T
... Ind2  # subpop. label: 6  (internally 1)
... C
... T
... END GENOTYPES
... ''')  
>>> outputfile = StringIO()
>>> convert2fasta(phasefile)
>Ind1_all1
T
>Ind1_all2
T
>Ind2_all1
C
>Ind2_all2
T
<BLANKLINE>


##### fastPhase file with a chromosome longer than the other ones (should give an error)
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T
... T C
... Ind2  # subpop. label: 6  (internally 1)
... C
... T
... END GENOTYPES
... ''')  
>>> convert2fasta(phasefile)
Traceback (most recent call last):
    ...
InvalidInputFile: Two chromosomes with different length
>>> 


##### fastPhase file with blank line
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T T T T T G A A A C C A A A G A C G C T G C G T C A G C C T G C A A T C T G
... 
... T T T T T G C C C C C A A A A G C G C G T C G T C A G T C T A A G A C C T A
... Ind2  # subpop. label: 6  (internally 1)
... C T T T T G C C C T C A A A A G T G C T G T G C C A G T C T A C G G C C T G
... T T T T T G A A A C C A A A G A C G C T T C G T C A G T A T A C G A T C T A
... END GENOTYPES
... ''')  

>>> convert2fasta(phasefile)
>Ind1_all1
TTTTTGAAACCAAAGACGCTGCGTCAGCCTGCAATCTG
>Ind1_all2
TTTTTGCCCCCAAAAGCGCGTCGTCAGTCTAAGACCTA
>Ind2_all1
CTTTTGCCCTCAAAAGTGCTGTGCCAGTCTACGGCCTG
>Ind2_all2
TTTTTGAAACCAAAGACGCTTCGTCAGTATACGATCTA
<BLANKLINE>

##### fastPhase file with MISSING SEQUENCE
>>> phasefile = StringIO('''
... ********************************************
... *                                          *
... *      Output from fastPHASE 1.2.3         *
... *      Code by P Scheet                    *
... *                                          *
... ********************************************
... BEGIN COMMAND_LINE
... /usr/bin/fastPHASE -T10 -K20 -p -u sample.popinfo -o sample sample.inp 
... END COMMAND_LINE
... 
... BEGIN COMMAND_EXPLAIN
...  K no. clusters (chosen or supplied): 20
...  S seed for random numbers (chosen or supplied): 1224252023
... END COMMAND_EXPLAIN
... 
... BEGIN DESCRIBE_TASKS
... minimize switch error
... END DESCRIBE_TASKS
... 
... BEGIN GENOTYPES
... Ind1  # subpop. label: 6  (internally 1)
... T T T T T G A A A C C A A A G A C G C T G C G T C A G C C T G C A A T C T G
... Ind2  # subpop. label: 6  (internally 1)
... C T T T T G C C C T C A A A A G T G C T G T G C C A G T C T A C G G C C T G
... T T T T T G A A A C C A A A G A C G C T T C G T C A G T A T A C G A T C T A
... END GENOTYPES
... ''')  

>>> convert2fasta(phasefile)
Traceback (most recent call last):
    ...
InvalidInputFile: Missing sequence in input file

"""
def main():
    pass

if __name__ == '__main__':
    import doctest
    doctest.testmod()