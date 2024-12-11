#!/usr/bin/env python3


import signal # for handling ctrl-c
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

from io_package import read_cell_and_pos_qe
from io_package import read_cell_and_pos_auto


def dQ(pre_gs,pre_es):
    """
    print the dQ between excited state and gs
    """
    (vecR, list_pos_f), package = read_cell_and_pos_auto(pre_gs)
    (vecR, list_pos_i), package = read_cell_and_pos_auto(pre_es)
    print("TODO")

def handler(signal_received, frame):
    # Handle any cleanup here
    sys.stderr.write('\nSIGINT or CTRL-C detected. Exiting gracefully\n')
    sys.exit(1)


def getPos(prefix):
    _, list_pos = read_cell_and_pos_qe(prefix)
    pos = np.zeros( (len(list_pos), 3) )
    for i, lp in enumerate(list_pos):
        pos[i] = lp['pos']
    # TODO use vecR to convert to angstrom
    return pos


def plotDensity(x, y, Z, color='Reds'):
    """
    3D plot of z(x,y)
    """
    # generate x,y mesh
    X, Y = np.meshgrid(x, y)

    # 
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # ax.contour3D(X, Y, abs(z), 50, cmap='binary')
    ax.plot_surface(X, Y, Z, cmap=color, linewidth=0, antialiased=False)

    ax.set_xlabel(r'x (\AA)')
    ax.set_ylabel(r'y (\AA)')
    ax.set_zlabel(r'$\Delta Q$')
    # ax.set_title('Direct')
    plt.savefig('density.png', dpi=300)

    return None


def main():
    """
    main program
    """
    # call handler if SIGINT or ctrl-c received
    signal.signal(signal.SIGINT, handler)

    pre_gs = "relax-gs/relax"
    pre_es = "relax-cdftup1/relax"

    pos_gs = getPos(pre_gs)
    pos_es = getPos(pre_es)

    deltaR = np.array( [ np.linalg.norm(g - e) for g, e in zip(pos_gs, pos_es) ] )

    # print(pos_gs[:,0])
    # print(deltaR)


    plotDensity(pos_gs[:,0], pos_gs[:,1], deltaR)


    
    return None


if __name__ == "__main__":
    main()

