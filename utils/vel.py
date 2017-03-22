import numpy as np
import matplotlib.pyplot as plt

empty_lines = 0
blocks = []
blocks_t = []
input_file = "output.dat"
for line in open(input_file):
	# Check for empty/commented lines
    if not line or line.startswith('#',0 , 1 ) \
                                        or line.startswith('\n',0 , 1) :
        # If 1st one: new block
        if empty_lines == 0:
            blocks.append([])
            empty_lines += 1
        # Non empty line: add line in current(last) block
    else:
        empty_lines = 0
        line = line.strip()
        blocks[-1].append(line.split())
    
    if line.startswith('#Time = '):
		line = line.split()
		blocks_t.append(line[2])
del blocks[-1]

blocks_t = np.asarray(blocks_t, np.float)

kB = 1. # Boltzman's constant
entropy = np.array([])
for block in blocks:
	blk = np.asarray(block, np.float64)

	vel = np.linalg.norm(blk[:, 4:7], axis = 1)

	bins = np.arange(np.floor(vel.min()),np.ceil(vel.max()))
	hist, bin_edges = np.histogram (vel, bins, density = True)

	entropy = np.append(entropy, -kB * np.dot(hist[np.nonzero(hist)],\
										np.log(hist[np.nonzero(hist)])))

plt.figure(0)
plt.plot(blocks_t, entropy)
plt.xlabel(r'Time, $\mathit{t}$')
plt.ylabel(r'Entropy, $\Delta \mathit{S}$')
plt.grid()
plt.savefig('S_t.pdf', format='pdf')

plt.figure(1)
plt.hist(vel, bins, normed = True)
plt.title("")
plt.xlabel(r'Velocity, $\mathit{v}$')
plt.ylabel(r'Normalized number of particles, $\Delta \mathit{N/N}$')
plt.grid()
plt.savefig('v_hist.pdf', format='pdf')

