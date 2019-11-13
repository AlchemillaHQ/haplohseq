import argparse
import sys
import random

#####################################################################
# This script takes in 2 input files and generates a pairwise LD map
# that can be used in pairwise phasing.
#
# Authors: Anthony San Lucas and Paul Scheet
#####################################################################
        
##### example snippet: snp pos file (can be a vcf file) #############
# chr22    16050408   
# chr22    16050612   
# chr22    16050678  
# chr22    16050984  
##################################################################### 

##### example snippet: hap file #####################################
# HG00096 HAPLO1 TCCCCTGACTATCGCTACAGCACTATGAGGGGAGAC
# HG00096 HAPLO2 TGCCCTGACTATCACTGCGGGGCTATGTGGGGAGAA
# HG00097 HAPLO1 TCCCCCCCCTGTCGCTACAGGGCTATGAGGGGAGAA
# HG00097 HAPLO2 CGTCATGACCATTGCTGCGGGGCCAAGATGGTAGGC
#####################################################################

class MapMaker:
    '''
    This class makes an LD Map given a set of haplotypes and a file specifying
    coordinates of the haplotype marker genomic positions.
    '''
    def __init__(self, markerFilename, haplotypesFilename):
        print "reading haplotype base coordinates ...."
        self.coords, self.coordIndexes, self.refs, self.alts = self.readMarkers(markerFilename)
#         self.haplotypes = self.readHaplotypes(haplotypesFilename)
        print "reading haplotype panel ...."
        self.markerCalls = self.readHaplotypes(haplotypesFilename, transpose = True)
    
    # reads from a vcf file of defined positions representing the haplotype data    
    def readMarkers(self, markerFilename):
        markerFile = open(markerFilename, 'r')
        coords = []
        refs = []
        alts = []
        coordIndexes = {}
        index = 0
        for line in markerFile:
            if line.strip().startswith("#"):
                continue
            tokens = line.strip().split('\t')
            coord = ":".join([tokens[0],tokens[1]])
            coords.append(coord)
            coordIndexes[coord] = index
            refs.append(tokens[3])
            alts.append(tokens[4])
            index += 1
        return coords, coordIndexes, refs, alts
    
    # reads from a file of predetermined haplotypes
    def readHaplotypes(self, haplotypesFilename, transpose):
        haplotypesFile = open(haplotypesFilename, 'r')

        if transpose is False:
            haps = []
            numHaps = 0
            for line in haplotypesFile:
                haps.append(list(line.strip().split()[2]))
    
            return haps
        else:
            markers = []
            firstHaplotype = haplotypesFile.readline().strip().split()[2]
            numMarkers = len(firstHaplotype)
            for index in range(0,numMarkers):
                markers.append([])
                markers[index].append(firstHaplotype[index])
                
            for line in haplotypesFile:
                haplotype = line.strip().split()[2]
                for index in range(0,numMarkers):
                    markers[index].append(haplotype[index])
            return markers
    
    # A and B are character vectors for site 1 and 2
    def dCalc(self, A, B, ref1, ref2, pR1, pR2, pA1, pA2):
        numRR = 0
        numRA = 0
        numHaps = len(A)
        
        for i in range(0, len(A)):
            if A[i] == ref1 and B[i] == ref2:
                numRR += 1
            elif A[i] == ref1 and B[i] != ref2:
                numRA += 1
        
        pRR = 1.0 * numRR / numHaps
        pRA = 1.0 * numRA / numHaps     
        
        dRR = pRR - pR1 * pR2        
        dRA = pRA - pR1 * pA2
        
        return dRR, dRA
        
    def genMap(self, pairDepth, log):
        # march through markers
        # collect ref and alt alleles and their frequencies
        # for every site, determine subsequent paired alleles
        
        polymorphicIndexes = []
        refAlleleCounts = []
        altAlleleCounts = []
        refAllelePairingMap = []
        numNonpolymorphicSites = 0
        
        polymorphicCoords = []
        polymorphicRefs = []
        polymorphicAlts = []
        
        # 1) COLLECT MARGINALS AND IDENTIFY SITES THAT ARE POLYMORPHIC
        # If a site has only one allele in the haplotype panel it is not informative, so it won't be used
        # in LD calculations.  We may want to consider setting some minimum threshold for a minor allele frequency.
        print "collecting marginal allele counts and identifying polymorphic sites from haplotype panel ...."
        for index in range(0,len(self.markerCalls)):
            ref = self.refs[index]
            alt = self.alts[index]
            
            refAlleleCount = 0
            altAlleleCount = 0
            
            markerAlleles = self.markerCalls[index]
            
            refAllelePairingFrequencies = []
            
            for allele in markerAlleles:
                if allele == ref:
                    refAlleleCount += 1
                else:
                    altAlleleCount += 1
                    
            if refAlleleCount == 0 or altAlleleCount == 0:
                numNonpolymorphicSites += 1
                continue
            else:
                polymorphicIndexes.append(index)
                polymorphicCoords.append(self.coords[index])
                polymorphicRefs.append(ref)
                polymorphicAlts.append(alt)
                
                refAlleleCounts.append(refAlleleCount)
                altAlleleCounts.append(altAlleleCount)
        
        log.write("Num polymorphic sites: " + str(len(polymorphicIndexes)) + "\n")
        log.write("Num non-polymorphic sites: " + str(numNonpolymorphicSites) + "\n")
               
        # 2) CALCULATE D FOR POLYMORPHIC SITES COMPARED WITH N SITES TO THE RIGHT THAT ARE POLYMORPHIC
        print "calculating D between polymorphic sites and their neighbors ...."
        dVals = []
        for i in range(0, len(polymorphicIndexes)):   # the last coordinate won't have any partners
            index_i = polymorphicIndexes[i]
            anchorCoord = polymorphicCoords[i]
            dVals.append([])
#             sys.stdout.write(anchorCoord + '\t')
            for j in range(i+1, min(len(polymorphicIndexes), i+1 + pairDepth)):
                index_j = polymorphicIndexes[j]
                partnerCoord = polymorphicCoords[j]
                refCount1 = refAlleleCounts[i]
                refCount2 = refAlleleCounts[j]
                altCount1 = altAlleleCounts[i]
                altCount2 = altAlleleCounts[j]
                
                ref1 = polymorphicRefs[i]
                ref2 = polymorphicRefs[j]
                alt2 = polymorphicAlts[j]
                
                pR1 = 1.0 * refCount1 / (refCount1 + altCount1)
                pR2 = 1.0 * refCount2 / (refCount2 + altCount2)
                pA1 = 1.0 * altCount1 / (refCount1 + altCount1)
                pA2 = 1.0 * altCount2 / (refCount2 + altCount2)
                
                dRR, dRA = self.dCalc(self.markerCalls[index_i], self.markerCalls[index_j], ref1, ref2, pR1, pR2, pA1, pA2)
                
#                 print "\t".join([str(dRR),str(dRA)])
                                
                if dRR >= dRA:
                    dVals[i].append(ref2)
                else:
                    dVals[i].append(alt2)
        
        return LDMap(polymorphicCoords, polymorphicRefs, polymorphicAlts, dVals, pairDepth)       
        
class LDMap:
    '''
    This class represents an LD Map
    '''
    def __init__(self, coords, chrCoords, refs, alts, dVals, pairDepth):
        self.coords = coords
        self.chrCoords = chrCoords  # add 1Mb bins to each chromosome
        self.coordIndexes = self.getIndexMap(coords)
        self.refs = refs
        self.alts = alts
        self.dVals = dVals
        self.pairDepth = pairDepth
    
    def getIndexMap(self, vec):
        elementIndexes = {}
        index = 0
        for element in vec:
            elementIndexes[element] = index
            index += 1
        return elementIndexes
    
    @staticmethod
    def fromFile(filename, pairDepth):
        chrCoords = {}
        coords = []
        refs = []
        alts = []
        dVals = []
        
        ldMapFile = open(filename, 'r')
        ldMapFile.readline()
        for line in ldMapFile:
            tokens = line.strip().split('\t')
            coordTokens = tokens[0].split(":")
            chr = coordTokens[0]
            coord = coordTokens[1]
            if chrCoords.get(chr) is None:
                chrCoords[chr] = []
            coords.append(tokens[0])
            chrCoords[chr].append(coord)
            refs.append(tokens[1])
            alts.append(tokens[2])
            if len(tokens) == 4:
                if pairDepth == 0:
                    pairDepth = len(tokens[3])
                if len(tokens[3]) <= pairDepth:
                    dVals.append(list(tokens[3]))
                else:
                    dVals.append(list(tokens[3])[0:pairDepth])

            else:
                dVals.append([])
        
        ldMap = LDMap(coords, chrCoords, refs, alts, dVals, pairDepth)
        return ldMap
        
    def save(self, filename):
        print "saving ldmap to " + filename
        ldMapFile = open(filename, 'w')
        ldMapFile.write("\t".join(["COORD","REF","ALT","PAIRED_ALLELES"]) + "\n")
        for i in range(0,len(self.coords)):
            ldMapFile.write("\t".join([self.coords[i],self.refs[i],self.alts[i],"".join(self.dVals[i])]) + "\n")
    
#     A = ['A','A','A','A','A','B','B','B']
#     B = ['B','A','A','A','B','A','A','A']
#     
#     mm = MapMaker(markerFilename, haplotypesFilename, mapFilename)
#     dRR, dRA = mm.dCalc(A,B,'A','A',5.0/8,6.0/8,3.0/8,2.0/8)

def generate_ldmap(log, pairDepth, markerFilename, haplotypesFilename, mapFilename):
    mm = MapMaker(markerFilename, haplotypesFilename)
    ldMap = mm.genMap(pairDepth, log)
    ldMap.save(mapFilename)   

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='''This script takes in 2 input files and generates a pairwise LD map
                                                    that can be used in pairwise phasing.''')
    parser.add_argument('-mp', '--marker_positions', 
                        help='''file of marker positions comprising the alleles of 
                                the haplotypes in the haplotype panel''', required=True)
    parser.add_argument('-hp', '--haplotype_panel', 
                        help='''file of haplotypes to use for LD calculations''', required=True)    
    parser.add_argument('-pd', '--pair_depth', 
                        help='''number of neighboring het sites to calculate LD values for (default=30)''', 
                        type = int,
                        default = 30)  
    parser.add_argument('-o', '--out_prefix', 
                        help='''output prefix for generated files''', required=True)   
    
    args = parser.parse_args()    
    log_filename = args.out_prefix + ".log"
    map_filename = args.out_prefix + ".ldmap"
    
    log = open(log_filename, "w")
    generate_ldmap(log, args.pair_depth, args.marker_positions, args.haplotype_panel, map_filename)



