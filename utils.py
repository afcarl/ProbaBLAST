import random
import numpy as np 

#----- Default values -----#
mismatch = -1
match = 2
gap = -2

#----- Methods ------#
def get_genome_and_probs():
	with open('Data/genome.txt', 'r') as handle:
		genome = list(handle.read().strip())
	with open('Data/probabilities.txt', 'r') as handle:
		probs = [float(f) for f in handle.read().strip().split(' ')]
	return genome, probs

def create_test_set(size=1000, length=25, sub_rate=0.003, ins_rate=0.003, del_rate=0.003):
	'''
	Creates a test test of size sequences of length length. These are 
	created by sampling repeatedly from the given genome.
	'''
	genome, probs = get_genome_and_probs()
	output = []
	ys = []
	for _ in xrange(0, size):
		ix = random.randrange(0, len(genome)-length)
		x = []
		mutations = 0
		for i in xrange(ix, ix+length):
			ns = ['A', 'C', 'T', 'G']
			d = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
			n, p = genome[i], probs[i]
			for t in ns:
				if t == n:
					d[t] = p
				else:
					d[t] = (1.-p)/3.
			x.append(np.random.choice(d.keys(), p=d.values()))
			dice = np.random.uniform()
			mutations += 1
			if dice < sub_rate:
				x[-1] = np.random.choice(ns)
			elif dice < ins_rate + sub_rate:
				x.append(np.random.choice(ns))
			elif dice < del_rate + sub_rate + ins_rate:
				del x[-1]
			else:
				mutations -= 1
			print mutations
		output.append(''.join(x)+','+str(ix))
	with open('test_%d.txt' % length, 'w') as handle:
		handle.write('\n'.join(output))

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




