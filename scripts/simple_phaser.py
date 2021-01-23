import argparse
import random
import sys

from ldmap import LDMap

#####################################################################
# This script takes in 2 input files: an ldmap and a vcf file.
# It generates candidate haplotypes for heterozygous sites in the
# vcf file that intersect sites in the ldmap.
#
# Authors: Anthony San Lucas and Paul Scheet
#####################################################################

class SimplePhaser:
    def __init__(self, ldMap):
        self.ldMap = ldMap
        
    def readPhasedHaps(self, machPhasedFilename):
        haps = []
        phasedFile = open(machPhasedFilename, 'r')
        #print machPhasedFilename
        hap1 = list(phasedFile.readline().strip().split()[2])
        hap2 = list(phasedFile.readline().strip().split()[2])
        haps.append(hap1)
        haps.append(hap2)
        return haps
    
    def readCoords(self, coordFilename):
        coords = []
        coordFile = open(coordFilename, 'r')
        for line in coordFile:
            if line.startswith("#"):
                continue
            tokens = line.strip().split("\t")
            if len(tokens) == 1:
                coords.append(tokens[0])
            else:
                coords.append(":".join([tokens[0],tokens[1]]))
        return coords
    
    def filterHaps(self, phasedHaps, coords, targetCoords):
        filteredHaps = []
        filteredHaps.append([])
        filteredHaps.append([])
        
        for index in range(0, len(coords)):
            if coords[index] in targetCoords:
                filteredHaps[0].append(phasedHaps[0][index])
                filteredHaps[1].append(phasedHaps[1][index])
        return filteredHaps
    
    def extractHetsFromVcf(self, vcfFilename, target_chrom=None):
        hets = []
        coords = []
        informativeVcfIndexes = []
        
        vcfFile = open(vcfFilename, 'r')
        
        i = 0
        previousCoord = None
        for line in vcfFile:
            if line.startswith("#"):
                continue
            tokens = line.strip().split('\t')
            chr = tokens[0]
            if target_chrom is not None and chr != target_chrom:
                continue
            if not chr.startswith("chr"):
                chr = "chr" + chr
            pos = tokens[1]
            ref = tokens[3]
            alt = tokens[4]
            genotype = tokens[9].split(':')[0]
            coord = ":".join([chr,pos])
            if previousCoord is None:
                previousCoord = coord
            else:   # skip third alleles
                if previousCoord == coord:
                    continue
                else:
                    previousCoord = coord
            if genotype == "0/0" or genotype == "1/1" or genotype == "./.":
                continue
            else:
                hets.append([ref,alt])
                coords.append(coord)
                informativeVcfIndexes.append(i)
            i += 1
                
        return hets, coords
    
    def getInformative(self, hets, coords):
        informativeHets = []
        informativeCoords = []
        informativeIndexes = []     # indexes into ldMap
        
        for i in range(0, len(coords)):
            tokens = coords[i].split(":")
            chr = tokens[0]
            pos = tokens[1]
            if self.ldMap.chrCoords.get(chr) is None:
                continue
            if pos in self.ldMap.chrCoords[chr]:
                informativeCoords.append(coords[i])
                informativeHets.append(hets[i])
                informativeIndexes.append(self.ldMap.coordIndexes[coords[i]])
        return informativeHets, informativeCoords, informativeIndexes
    
    def getRefPairedAllele(self, anchorIndex, pairedIndex):
#         print("\t".join([str(anchorIndex), str(pairedIndex - anchorIndex  - 1), str(len(self.ldMap.dVals)),str(len(self.ldMap.dVals[anchorIndex]))]) # prints intermarker counts)
        if (pairedIndex - anchorIndex - 1) >= len(self.ldMap.dVals[anchorIndex]):
            pRef = random.uniform(0,1)
            if pRef > 0.5:
                return self.ldMap.refs[pairedIndex]
            else:
                return self.ldMap.alts[pairedIndex]
            
        return self.ldMap.dVals[anchorIndex][pairedIndex - anchorIndex - 1]
        
    
    def phaseVcf(self, vcfFilename, target_chrom=None):
        hets, coords = self.extractHetsFromVcf(vcfFilename, target_chrom)
        phasedAlleles, informativeCoords = self.phase(hets, coords)
        return phasedAlleles, informativeCoords
            
    def phase(self, hets, coords):
        
        informativeHets, informativeCoords, informativeIndexes = self.getInformative(hets, coords)
        
        #print "informativeHets " + str(len(informativeHets))
        #print "informativeCoords " + str(len(informativeCoords))
        #print "informativeIndexes " + str(len(informativeIndexes))
        
        phasedAlleles = []
        phasedAlleles.append([])
        phasedAlleles.append([])
        
        # initially, the reference haplotype (the anchorHap) is 0
        phasedAlleles[0].append(self.ldMap.refs[informativeIndexes[0]])
        phasedAlleles[1].append(self.ldMap.alts[informativeIndexes[0]])
        previousCoord = coords[0]
        anchorHap = 0
        altHap = 1

        for i in range(1,len(informativeIndexes)):
            anchorIndex = informativeIndexes[i-1]
            pairedIndex = informativeIndexes[i]
            
            pairedAllele = self.getRefPairedAllele(anchorIndex, pairedIndex)
            
            ref = self.ldMap.refs[informativeIndexes[i]]
            alt = self.ldMap.alts[informativeIndexes[i]]
            
            if pairedAllele == ref:
                phasedAlleles[anchorHap].append(ref)
                phasedAlleles[altHap].append(alt)
            
            else: # assume pairedAllele == alt:
                phasedAlleles[anchorHap].append(alt)
                phasedAlleles[altHap].append(ref)
                
                # swap haps
                temp = anchorHap
                anchorHap = altHap
                altHap = temp 
            
            if pairedAllele not in [ref,alt]:
                print("pairedAllele not ref or alt " + ref + " " + alt)
                print(pairedAllele)
         
        return phasedAlleles, informativeCoords
            
        
    def phaseConcordance(self, phasedHaps, candidateHap):
        concordances = []
        previousMatch = 0

        for i in range(0, len(candidateHap)):
            if candidateHap[i].upper() == phasedHaps[0][i].upper():
                currentMatch = 0
            elif candidateHap[i].upper() == phasedHaps[1][i].upper():
                currentMatch = 1
            else:
                print("ERROR: " + candidateHap[i] + " does not match " + phasedHaps[0][i] + " or " + phasedHaps[1][i] + "\t" + self.ldMap.coords[i] + "\t" + self.ldMap.refs[i] + "\t" + self.ldMap.alts[i])
#                 # randomly assign match
#                 p0 = random.uniform(0,1)
#                 if p0 > 0.5:
#                     currentMatch = 0
#                 else:
#                     currentMatch = 1
                continue
            
            if i != 0 and currentMatch == previousMatch:
                concordances.append(1)
            else:
                concordances.append(0)
            previousMatch = currentMatch
        concordance = float(sum(concordances))/len(concordances) if len(concordances) > 0 else float('nan')
        #print "num concordance points: " + str(len(concordances))
        return concordance

    def writeHaps(self, haps, outfile):
        '''
        This method writes out pairwise-phased haplotypes.
        '''
        hapFile = open(outfile, "w")
        simpleName = outfile.split("/")[-1]
        hapFile.write(" ".join([simpleName,"HAP0", "".join(haps[0])]) + "\n")
        hapFile.write(" ".join([simpleName,"HAP1", "".join(haps[1])]) + "\n")
        hapFile.close() 
        
    def filterVcf(self, vcfFilename, phasedCoords):
        '''
        Given a vcf file and a set of coordinates, this method will return the
        vcf lines that intersect the given coordinates.
        '''
        vcfOutLines = []
        vcfFile = open(vcfFilename, 'r')
        numGT = 0
        previousCoord = None        # This logic to remove third alleles shouldn't be needed.
        
        for line in vcfFile:
            if line.startswith("#"):
                vcfOutLines.append(line.rstrip("\n"))
            else:
                tokens = line.rstrip("\n").split("\t")
                chr = tokens[0]
                pos = tokens[1]
                if not chr.startswith("chr"):
                    chr = "chr" + chr
                coord = ":".join([chr,pos])
                if previousCoord is None:
                    previousCoord = coord
                else:
                    if previousCoord == coord:
                        continue    # only allow 1 variant per site
                    previousCoord = coord
                if coord in phasedCoords:
                    vcfOutLines.append(line.rstrip("\n"))
                    numGT += 1
        return vcfOutLines, numGT
        
    def writeVec(self, vec, delim, outFilename):
        outFile = open(outFilename, "w")
        outFile.write(delim.join(vec) + "\n")
        outFile.close()

    def writePos(self, pos, outfile):
        '''
        This method writes out positions of the haps.
        '''
        posFile = open(outfile, "w")
        posFile.write("\n".join(pos))
        posFile.close()        

def run(log, vcfFilename, mapFilename, pairDepth, outPrefix, target_chrom=None):

    log.write("loading LD map " + mapFilename + " with pair depth " + str(pairDepth) + "\n")
    print("loading LD map " + mapFilename + " with pair depth " + str(pairDepth))
    ldMap = LDMap.fromFile(mapFilename, pairDepth)
    log.write("loaded ldmap " + mapFilename + "\n")
    
    pp = SimplePhaser(ldMap)
    print("phasing " + vcfFilename)
    log.write("phasing " + vcfFilename + "\n")
    phasedHaps, phasedCoords = pp.phaseVcf(vcfFilename, target_chrom)
    log.write("num phased coords (" + vcfFilename + "): " + str(len(phasedCoords)) + "\n")
    
    filteredVcf, numGT = pp.filterVcf(vcfFilename, phasedCoords)
    log.write("num vcf genotypes (" + vcfFilename + "): " + str(numGT) + "\n")
    pp.writeVec(filteredVcf, "\n", outPrefix + ".hap.vcf")
    pp.writeHaps(phasedHaps, outPrefix + ".hap")
    pp.writePos(phasedCoords, outPrefix + ".pos")
    print("generated:\n\t" + outPrefix + ".hap\n\t" + outPrefix + ".pos\n\t" + outPrefix + ".hap.vcf")
#     print len(phasedHaps[0])
#     print len(phasedHaps[1])
#     print len(phasedCoords)
#     print numGT

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='''This script takes in 2 input files: an ldmap and a vcf file.
                                                    It generates candidate haplotypes for heterozygous sites in the
                                                    vcf file that intersect sites in the ldmap.''')
    parser.add_argument('-ld', '--ldmap', 
                        help='''ldmap filename''', required=True)
    parser.add_argument('-v', '--vcf', 
                        help='''vcf file of markers for phasing''', required=True)    
    parser.add_argument('-pd', '--pair_depth', 
                        help='''number of neighboring het sites to calculate LD values for (default=30)''', 
                        type = int,
                        default = 30)  
    parser.add_argument('-tc', '--target_chrom', 
                        help='''specify if want to limit to 1 chromosome''')  
    parser.add_argument('-o', '--out_prefix', 
                        help='''output prefix for generated files''', required=True)   
    
    args = parser.parse_args()  
    
    log = open(args.out_prefix + ".log", "w")
    run(log, args.vcf, args.ldmap, args.pair_depth, args.out_prefix, args.target_chrom)
    log.close()
          
