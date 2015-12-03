

#----- Default values -----#
mismatch = -1
match = 2
gap = -2

#----- Methods ------#

def score(S,Ps,T):
	"""
	Calculates alignment score of sequence S and T with respect to probabilities for S 
	Matches: +2
	Mismatches: -1
	Gaps: -2

	Parameters
	----------
	S:	array-like or string, one character at each position
	Ps:	array-like, probability in [0,1] that letter at position i was correctly read
	T:	array-like or string, one character at each position
	Returns
	-------
	score: float 
	"""

	assert(len(S)==len(T))
	score = 0
	
	for s,p,t in zip(S,Ps,T):
		if s == t:
			score += match*p 
		elif s=='-' or t=='-':
			score += gap 
		else:
			score += mismatch*p 

	return score



