import numpy as np
from PointinPolygon import inpolygon
import os
import glob
import h5py

# Generate Target Points matrix X,Y
nx = 300
ny = 400
x = np.linspace(-98.12,-97.735,nx)
y = np.linspace(36.59,37.00,ny)
X,Y = np.meshgrid(x,y)
X = X.reshape(-1)
Y = Y.reshape(-1)
n = len(X)
xv = [-98.11,-98.065,-98.065,-98.05,-97.95,-97.745,-97.745,-97.875,-98.07,-98.11,-98.11]
yv = [36.60,36.60,36.63,36.63,36.72,36.765,36.99,36.99,36.82,36.82,36.60]

# Set the minima sample times
min_number = 20

dsfiles = glob.glob('outputs_random/*h5')


## reading the Voronoi Partitions
files = []
for i in range(n):
	files.append([])
for j in range(len(dsfiles)):
	h5file = h5py.File(dsfiles[j],'r')
	xvyv = h5file['xvyv'][:]
	h5file.close()
	_in,_on = inpolygon(X, Y, xvyv[0,:], xvyv[1,:])
	indx = np.array([i for i,x in enumerate(_in) if x])
	for i in indx:
		files[i].append(j)

## reading the dispersion curves
curves = []
for fr in dsfiles:
	fr1 = os.path.join('curves/',fr.split('/')[-1].split('.')[0]+'.txt')
	curve = np.loadtxt(fr1)
	curves.append(curve)

outdir = 'data'
if not os.path.exists(outdir):
	os.mkdir(outdir)

# estimate the dispersion curves of each target points
for i in range(n):
	print(i)
	_in,_on = inpolygon(X[i],Y[i],xv,yv)
	if not _in:
		continue
	fid = open(os.path.join(outdir,'curve'+str(i)+'.txt'),'w')
	for modei in range(1):
		cs = []
		for j in files[i]:
			curve = curves[j]
			if len(curve[curve[:,2]==modei,2])>0:
				cs.append(curve[curve[:,2]==modei,:])
	
		minf = 10
		maxf = 0
		for csi in cs:
			if minf>np.min(csi[:,0]):
				minf = np.min(csi[:,0])

			if maxf <np.max(csi[:,0]):
				maxf = np.max(csi[:,0])
		nf = int((maxf-minf)/0.1) + 1
		for k in range(nf):
			c = []
			f = np.round(minf + k*0.1,2)
			for csi in cs:
				indx = np.where(np.abs(csi[:,0]-f)<0.01)[0]
				if len(indx) > 0:
					c.append(csi[indx,1])
			c = np.asarray(c)/1e3
			if len(c) > min_number:
				mean = np.round(np.mean(c),3)
				std = np.round(np.std(c),3)
				fid.write(str(f)+' '+str(mean)+' '+str(modei)+' '+str(std)+'\n')
	fid.close()