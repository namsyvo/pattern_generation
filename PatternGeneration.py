"""
Generating patterns with Wilcoxon rank sum test
"""

import Wilcoxonranksum


"""
Generating all set of n elements from a list (items)
"""

def combinations(items, n):
    if n == 0: yield []
    else:
        for i in range(len(items)):
            for cc in combinations(items[i+1:], n - 1):
                yield [items[i]] + cc


"""
Creating patterns
"""

def gen_pattern(alpha, *treatment_groups):

	CombList = [] 
	for y in combinations(range(1,(len(treatment_groups) + 1)), 2): CombList.append(tuple(y))

	treatdict = {}
	pattern = ''

	for m in range(len(treatment_groups)):
		treatdict[m + 1] = treatment_groups[m]

	for (a,b) in CombList:
		Pvalue, PvalueC = Wilcoxonranksum.RanksumComparison(treatdict[a], treatdict[b])
		#print treatdict[a], treatdict[b], Pvalue, PvalueC
		if Pvalue < alpha:
			pattern = pattern + str(2)
		elif PvalueC < alpha:
			pattern = pattern + str(0)
		else:
			pattern = pattern + str(1)
	return pattern
