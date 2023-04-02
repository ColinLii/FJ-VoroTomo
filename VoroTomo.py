from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import pandas as pd
from PointinPolygon import inpolygon
import os
from geopy.distance import great_circle
import ccfj
import h5py
from scipy import stats


def Pairs(sta):
    p = []
    nsta = len(sta)
    for ii in range(nsta):
        for jj in range(ii+1,nsta):
            p.append([sta[ii],sta[jj]])
    return p
def cal_indx(pair,nsta):
    indx = int(pair[0]*(2*nsta-pair[0]-1)/2+pair[1]-pair[0]-1)
    return indx

def singleVoro(indx0,outname,xvyv):
    global ncfs
    global f
    global c
    global nf
    global r
    global count
    global outdir
    filename = os.path.join(outdir,outname+'.h5')
    if os.path.exists(filename):
        os.remove(filename)
    
    subpairs = Pairs(indx0)
    indx1 = [cal_indx(pair,nsta) for pair in subpairs]
    ncfsi = ncfs[indx1,:]
    counti = count[indx1]
    ri = r[indx1]
    indx = np.argsort(ri)
    ri = ri[indx]
    ncfsi = ncfsi[indx,:]
    counti = counti[indx]
    ncfsi = ncfsi[counti!=0]
    ri = ri[counti!=0]
    counti = counti[counti!=0]
    for i in range(len(ri)):
        ncfsi[i,:] = ncfsi[i,:]/counti[i]

    
    ds = ccfj.fj_noise(np.real(ncfsi),ri*1e3,c,f,itype=1,func=1)

    h5file=h5py.File(filename,'w')
    h5file.create_dataset('ds',data=ds)
    h5file.create_dataset('c',data=c)
    h5file.create_dataset('f',data=f)
    h5file.create_dataset('indx',data=indx0)
    h5file.create_dataset('xvyv',data=xvyv)
    h5file.close()


if __name__ == '__main__':
    # input Dir
    Dir = './'
    # output Dir
    outdir = 'outputs_random'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # excel file saving the location information
    stainfo =  pd.read_excel(os.path.join(Dir,'sta.xlsx'))
    # np.array([longitude, latitude])
    staloc =  np.asarray([stainfo.iloc[:,6],stainfo.iloc[:,5]]).T
    # The KDE of the stations' locations
    kernel = stats.gaussian_kde(staloc.T)
    nsta = len(stainfo.iloc[:,0])
    # read the NCFs
    ncffile = h5py.File(os.path.join(Dir,'ncfs_g1.h5'),'r')
    ncfs = ncffile['ncfs'][:]
    r = ncffile['r'][:]
    count = ncffile['count'][:]
    f = ncffile['f'][:]
    nf = len(f)

    # Phase velocity range
    c = np.linspace(800,4000,1000)

    # VoroTomo Partition times
    n = 100
    # number of Voro cells  range
    kmin = 20
    kmax = 40

    # Set a larger range of the stations to facilitate the calculation of Voronoi cells
    # the points0 are not used in the inner Voronoi cells generating
    Radius = 10 
    x0 = np.mean(staloc[:,0])
    y0 = np.mean(staloc[:,1])
    minx = min(staloc[:,0]);maxx = max(staloc[:,0])
    miny = min(staloc[:,1]);maxy = max(staloc[:,1])
    points0 = []
    for i in range(36):
        points0.append([x0+np.cos(i/18*np.pi),y0+np.sin(i/18*np.pi)])
        
    points0=np.asarray(points0)

    # Starting Voronoi Partitioning
    for ii in range(n):
        # random the number of Voronoi cells
        k = np.random.randint(kmin,kmax)
        # Using the KDE PDF to generate Voronoi cells
        points = kernel.resample(k).T 
        # Using the Uniform distribution to generate Voronoi cells
        '''
        points = np.random.rand(k,2)
        points[:,0] = points[:,0]*(maxx-minx) + minx
        points[:,1] = points[:,1]*(maxy-miny) + miny
        '''

        # Voronoi Partition
        points = np.concatenate((points,points0))
        vor = Voronoi(points)
        areas = []
        for j in range(k):
            xv = vor.vertices[vor.regions[vor.point_region[j]],0]
            yv = vor.vertices[vor.regions[vor.point_region[j]],1]
            _in,_on = inpolygon(staloc[:,0],staloc[:,1],xv,yv)
            indx = np.array([i for i,x in enumerate(_in) if x])
            if len(indx)>10:
                areas.append([indx,np.array([xv,yv])])
        vor.close()
        for j in range(len(areas)):

            print(ii,j,len(areas[j][0]))
            # calculate the F-J spectrum
            singleVoro(areas[j][0],'vor'+str(ii)+'_'+str(j),areas[j][1])