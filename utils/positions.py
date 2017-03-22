#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# Save figures to single pdf
#from matplotlib.backends.backend_pdf import PdfPages
import math

# Read data blocks
empty_lines = 0
blocks = []
blocks_t = []
input_file = sys.argv[1]
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

# Read boundaries blocks
empty_lines = 0
boundaries = []
boundaries_name = []
border_file = sys.argv[2]
for line in open(border_file):
	# Check for empty/commented lines
    if not line or line.startswith('#',0 , 1 ) \
                                        or line.startswith('\n',0 , 1) :
        # If 1st one: new block
        if empty_lines == 0:
            boundaries.append([])
            empty_lines += 1
        # Non empty line: add line in current(last) block
    else:
        empty_lines = 0
        line = line.strip()
        boundaries[-1].append(line.split())
    
    if line.startswith('#Boundary'):
		line = line.split()
		boundaries_name.append(line[1])
# Close boundary perimeter line
for boundary in boundaries:
	boundary.append(boundary[0])


#pdffile = PdfPages('positions.pdf')
fig = plt.figure()
# The following code is for setting equal scales of axes
ax = fig.add_subplot(111, projection='3d')
#pos = np.array([])
pos = np.asarray(blocks[0], np.float64)[:, 1:4]
max_range = np.array([pos[:, 0].max()-pos[:, 0].min(), \
					  pos[:, 1].max()-pos[:, 1].min(), \
					  pos[:, 2].max()-pos[:, 2].min()]).max() / 2.0

mid_x = (pos[:, 0].max()+pos[:, 0].min()) * 0.5
mid_y = (pos[:, 1].max()+pos[:, 1].min()) * 0.5
mid_z = (pos[:, 2].max()+pos[:, 2].min()) * 0.5 
# #

for cnt, block in enumerate(blocks):
	blk = np.asarray(block, np.float64)
	pos = blk[:, 1:4]
	col = blk[:, 0]

	# Plot title
	ax.set_title('Particles positions,\n t = ' \
					+ str(blocks_t[cnt]) + ' s')
	# Axes limits
	ax.set_xlim(mid_x - max_range, mid_x + max_range)
	ax.set_ylim(mid_y - max_range, mid_y + max_range)
	ax.set_zlim(mid_z - max_range, mid_z + max_range)
	# Axes labels
	ax.set_xlabel(r'$\mathit{x}$')
	ax.set_ylabel(r'$\mathit{y}$')
	ax.set_zlabel(r'$\mathit{z}$')
	
	# Plot boundaries
	for boundary in boundaries:
		bnd = np.asarray(boundary, np.float64)
		ax.plot(bnd[:, 0], bnd[:, 1], bnd[:, 2], \
		        marker = '+', markersize = 1.0, \
		        linestyle = '-', color = 'red', alpha = 0.2)
		
	# Plot particles
	ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c = col, \
				edgecolors = 'black', s = math.pi)
	#fig.tight_layout()
	fig.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0, hspace=0)

	
	print('Plotting block #' + str(cnt) +\
		  '\t t = ' + str(blocks_t[cnt]) + ' s')
	
	# Save plot
	fig.savefig(str(cnt).zfill(6) + '.jpg', format='jpg')
	# Save figure to pdf
	#pdffile.savefig(fig)
	# Clear axis
	ax.cla()
	#Stop plotting at block number
	#if (cnt == 1):	break
	#cnt+=1

#pdffile.close()

