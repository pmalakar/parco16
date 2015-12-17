import sys
import numpy as np

times = []

f = open(sys.argv[1])

for line in f.readlines():
		words = line.split()
		times.append(float(words[15]))
	
f.close()

#print 'Min ', np.min(times), ' Max ', np.max(times), ' Mean ', np.mean(times), ' Median ', np.median(times)
print np.min(times), np.max(times), np.mean(times), np.median(times)

