import numpy as np
from matplotlib.path import Path
from matplotlib import pyplot as plt

def inpolygon(xq, yq, xv, yv):
    vertices = np.vstack((xv, yv)).T
    path = Path(vertices)
    test_points = np.hstack([xq.reshape(xq.size, -1), yq.reshape(yq.size, -1)])
    _in = path.contains_points(test_points)
    _in_on = path.contains_points(test_points, radius=-1e-10)
    _on = _in ^ _in_on
    
    return _in_on, _on

def incircle(xq,yq,xc,yc,radius):
    dist = np.array(np.sqrt((xq-float(xc))**2+(yq-float(yc))**2))
    indx=np.where(dist<radius)
    indx = indx[0]
    return indx

if __name__ == '__main__':
    xv = np.array([-4, 0, 4, 0])
    yv = np.array([0, 4, 0, -4])
    X = np.array([0, 1, 3.5, 4, 5])
    Y = np.array([0, 1, 0,   0, 0])
    xv = np.append(xv,xv[0])
    yv = np.append(yv,yv[0])
    _in, _on = inpolygon(X, Y, xv, yv)
    print(_in, _on)
    plt.plot(xv,yv)
    plt.plot(X,Y,'bo')
    indx = np.array([i for i,x in enumerate(_in) if x])
    plt.plot(X[indx],Y[indx],'ro')
    plt.show()
