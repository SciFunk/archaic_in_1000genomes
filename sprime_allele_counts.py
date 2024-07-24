import numpy as np
import pandas as pd
from collections import defaultdict
import sys
import gzip

##dictionary to tie individual to population
ppl_pop_dicc ={}
with open('pop_locations.txt', "r") as f:
    poplabels = []
    for line in f:
        spline = line.split()
        poplabels.append(spline[1])
        ppl_pop_dicc[spline[0]] = spline[1]
        poparray = np.asarray(poplabels)
        #Gives each unique item in array: just one time, in alphabetic order
        popnames = np.unique(poparray)

##file containing all rsIDs in sPrime data
IDTable = pd.read_csv('rsIDallpops_insPrime.csv')
posList = IDTable[['ID','arc_allele']].values.tolist()

##list of population (abbreviations) in 1000 Genomes data
poplist = ['LWK', 'GWD', 'MSL', 'YRI', 'ESN', 'ACB', 'ASW', "PEL", "MXL", "CLM", "PUR", "JPT", "CDX", "KHV", "CHB", "CHS", "FIN", "CEU", "GBR", "IBS", "TSI", "GIH", "ITU", "STU", "PJL", "BEB"]
poptotals = defaultdict(int)

##open each chromosome file and count allele frequency in each population
with gzip.open(sys.argv[1], 'rt') as f:
    for _ in range(252):
        next(f)
    for line in f:
        spline = line.split()
        if spline[0] == "#CHROM":
            header=spline
            pop_index_dicc = {}
            for i in range(9,len(header)):
                pop = ppl_pop_dicc[header[i]]
                try:
                    pop_index_dicc[pop].append(i)
                except KeyError:
                    pop_index_dicc[pop] = [i]
            for pop in poplist:
                poptotals[pop] = float(len(pop_index_dicc[pop]))
                print(poptotals[pop])
        else:
            rsID = str(spline[2])
            ref = str(spline[3])
            alt = str(spline[4])
            for item in posList:
                if str(item[0]) == rsID:
                    print(str(item[0]),rsID)
                    arrayspline = np.asarray(spline)
                    if str(item[1]) == ref:
                        for pop in poplist:
                            poplocations = arrayspline[pop_index_dicc[pop]]
                            count01 = len(np.where(poplocations=='0|1')[0])
                            count10 = len(np.where(poplocations=='1|0')[0])
                            count11 = 2.0*(len(np.where(poplocations=='0|0')[0]))
                            allcounts = count01 + count10 + count11
                            alleleFreq = allcounts/(2*poptotals[pop])
                            item.append(alleleFreq)
                        print(item)
                    if str(item[1]) == alt:
                        for pop in poplist:
                            poplocations = arrayspline[pop_index_dicc[pop]]
                            count01 = len(np.where(poplocations=='0|1')[0])
                            count10 = len(np.where(poplocations=='1|0')[0])
                            count11 = 2.0*(len(np.where(poplocations=='1|1')[0]))
                            allcounts = count01 + count10 + count11
                            alleleFreq = allcounts/(2*poptotals[pop])
                            item.append(alleleFreq)
                        print(item)


posListdf = pd.DataFrame(posList)
match = re.search(r'ALL.(chr\d+)', sys.argv[1])
posListdf["chrom"] = match.group(1)
posListdf.to_csv(path_or_buf = r"{}.csv".format(match.group(1)), mode='a')
#posListdf.to_csv("C:\\Users\\SciFunk\\Downloads\\Working\\Biostats\\MatchFiles\\NATEASallelefreq.csv", mode='a')
