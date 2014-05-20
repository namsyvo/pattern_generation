"""
Calculation exact Pvalue for Wilcoxon rank sum test using recursive method
"""

import math

""" 
Function SortnRank takes a list and returns sorted list and ranked list containing ranks of elements 
preserving order of that of original list.
Input: List containing data
Output: Sorted List , Ranked List 
"""

def SortnRank(Array):
	
	Arraysort = Array[:]
	Arraysort.sort(lambda a,b: cmp(float(a), float(b))) 
	Arraysortc = Arraysort[:] 
	RankedArray = [0]*len(Arraysortc)
	
	for m in Arraysortc:
		RankedArray[Arraysortc.index(m)] = Array.index(m)
		Arraysortc[Arraysortc.index(m)] = -1
		Array[Array.index(m)] = -1
	
	return Arraysort,RankedArray


"""
Function RanknTieresolver takes List and resolves ties in ranks and substitutes them with average value.
This function uses Function SortnRank for sorting and ranking.
Input: List containing data.
Output: Ranked list with ties resolved.
"""

def RanknTieresolver(inputlist):
    
    n = len(inputlist)
    sortlist, ranked = SortnRank(inputlist)
    sumranks = 0 
    dupcount = 0
    
    newlist = [0]*n
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i==n-1 or sortlist[i] != sortlist[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newlist[ranked[j]] = averank
            sumranks = 0
            dupcount = 0
    return newlist


"""
Function RanksumComparison prepares input to be sent to RecursionPvalue program. This function is based on 
Wilcoxon Rank Test procedure but uses Recursive method for obtaining P-values.
Input: Two Lists to be compared (List1, List2).
Output: Pvalue for two cases Mean(List1) < Mean(List2) (Pvalue1), Mean(List2) < Mean(List1) (Pvalue2).
"""

def RanksumComparison(List1, List2):

	CommonList = []
	RankedList = []
	
	Numsample1 = len(List1)
	Numsample2 = len(List2)
	CommonList = List1 + List2
	RankedList = RanknTieresolver(CommonList)
	RanksumList1 = sum(RankedList[0:Numsample1])

	"""
	Initialize P (memory array)
	"""
	depth = max(Numsample1, Numsample2)
	size = 2*depth**2 + int(RanksumList1) + 5
	P = [[[-1.0 for i in range(size)] for j in range(Numsample2+1)] for k in range(Numsample1 + 1)]

	"""
	Take index from r
	"""
	def idx(r):
		return int(r)+2*depth**2+1

	"""
	Reference: M.S. Ross Simulation 3rd edition, 2002.
	This function implements Recursive equation given by P(n,m,r) = n/(n+m)*P(n-1,m,r-n-m)+ m/(n+m)*P(n,m-1,r).
	Base cases: P(1,0,k<=0) = 0, P(1,0,k>0) = 1, P(0,1,k<0) = 0, P(0,1,k>=0) = 1.
	Argument n, m are positive integer types and r can be integer or float.
	Input: Elements in Sample 1, Elements in Sample 2, Sum of ranks of elements in Sample 1.
	Output: Pvalue
	"""

	def getPvalues(n, m, r):

		if n == 0.0:
			if m == 1.0:
				if r < 0.0:
					return 0.0
				else:
					return 1.0
	
		if n == 1.0:
			if m == 0.0:
				if r > 0.0:
					return 1.0
				else:
					return 0.0	
			
		if n == 0.0:
			return getPvalues(n, m-1, r)
		if m == 0.0:
			return getPvalues(n-1, m, r-n-m)
	
		if P[n][m][idx(r)] != -1.0:
			return P[n][m][idx(r)]
		else:
			P[n][m][idx(r)] = (float(n)/float(n+m))*(getPvalues(n-1, m, r-n-m)) +  (float(m)/float(n+m))*(getPvalues(n, m-1, r))
			return P[n][m][idx(r)]

	Pvalue = getPvalues(Numsample1, Numsample2, RanksumList1)

	"""
	Re-initialize P
	"""
	P = [[[-1.0 for i in range(size)] for j in range(Numsample2 + 1)] for k in range(Numsample1 + 1)]
	PvalueC = 1 -(getPvalues(Numsample1, Numsample2, (RanksumList1 - 1)))
	
	return Pvalue, PvalueC
