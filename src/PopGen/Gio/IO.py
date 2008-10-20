#!/usr/bin/env python
# Read various different Population Genetics formats (ped, tped, fam, etc..)

from Individual import Individual
from Population import Population

def read_tped(tped_path):
    """
    Parse a TPED file and returns a dictionary of population objects.
    
    """
    print tped_path

def read_ped(ped_path):
    """
    Parse a PED file and returns a dictionary of population objects.
    Every line in the PED file corresponds to an individual. However, we will only count the number of hetherozygotes for every marker for every population.

    Example of PED file:
    # population id father mother char1 phenotype locus1_allele1 locus1_allele2 locus2_allele1 .....
    Mandenka HGDP00912 0 0 1 2 C C G A G G T T T T A A G G
    Mandenka HGDP01283 0 0 1 2 T C G A G G C T C C G G G A
    Yoruba HGDP00928 0 0 2 2 C C G G G G C T T T G A G G
    Yoruba HGDP00937 0 0 1 2 C C A A A G C C T C A A G A

    Example of output:
    populations = {Mandenka: [HGDP00912, HGDP01283], Yoruba: [HGDP00928, HGDP00937]}
    """
    print "parsing ped file" + ped_path
    print

    populations = {}

    for line in file(ped_path, 'r'):
        if not line.startswith('#'):
            # skip line if it is a comment (starts with '#'). Should be replaced with python 2.5 'with' syntax.
            
            ped_fields = line.split()

            population_name = ped_fields[0]
            if not populations.has_key(population_name):
                # if this is the first time an individual of this population is found, create a new empty entry in the populations dictionary.
                populations[population_name] = Population(population_name)
                
            # Create an Individual object from the current line.
            current_ind = Individual().from_ped_line(ped_fields)
            print current_ind
#            print ped_fields[0:]

    return populations

            
if __name__ == '__main__':
    ped_test_path = '../../test/data/ped/malaria_test_set.ped'
    pops = read_ped(ped_test_path)
#    for pop in pops:
#        print pops[pop].name

