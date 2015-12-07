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

def create_test_set(size=1000, length=25, sub_rate=0.01, ins_rate=0.01, del_rate=0.01):
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


def make_convex(fpr, tpr):
	init_len = len(fpr)
	fpr_new = []
	tpr_new = []
	
	for i in range(init_len):
		if i == 0 or i == init_len-1:
			fpr_new.append(fpr[i])
			tpr_new.append(tpr[i])
		else:
			x_prev = fpr[i-1]
			y_prev = tpr[i-1]
			x = fpr[i]
			y = tpr[i]
			x_next = fpr[i+1]
			y_next = tpr[i+1]
			if x_prev == x and y_prev < y:
				fpr_new.append(fpr[i])
				tpr_new.append(tpr[i])
				continue
			elif x == x_next and x_next < x:
				fpr_new.append(fpr[i])
				tpr_new.append(tpr[i])
				continue
			# params y = ax + b
			a = (y_next-y_prev)/(x_next-x_prev)
			b = y_next - a * x_next
			y_on_line = a*x + b
			if y_on_line < y:
				fpr_new.append(fpr[i])
				tpr_new.append(tpr[i])

	if init_len > len(fpr_new):
		fpr_new, tpr_new = make_convex(fpr_new,tpr_new)

	return fpr_new, tpr_new

def remove_duplicates(rates):

	rates.sort()
	return list(rates for rates,_ in itertools.groupby(rates))


