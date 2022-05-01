"""
Animated snapshots
"""

from dovis.core import *

import numpy as np, matplotlib.pyplot as plt, sys
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize

def arrows(args):
    """
    Movie of flying arrows
    """
    block = process_block(args.file)

    # parameters will set the scene
    P = next(block)
    if(len(P.shape) > 2): raise RuntimeError("Can't (yet) create movies when d>2")

    # Backdrop

    if len(P.shape) == 1:
        # Determine aspect ratio
        L = P.shape[0]
        aspect = max(1,L/20)
        fig, ax = plt.subplots(1,1, figsize=(6,6/aspect))
        ax.set_xlim(0, L)
        ax.set_ylim(0,1)
        ax.set_title(f't={P.initial}')
        ax.set_aspect(aspect)
        q = ax.quiver( np.arange(L)+0.5, 0.5, np.ones(L), np.zeros(L), pivot='mid', cmap=plt.get_cmap('Greys'),norm=Normalize(0,1))
        # Update function
        def update(ns):
            n, S = ns
            dx = P.dirs[S.directions.flatten()]
            dy = np.zeros_like(dx)
            if args.active:
                q.set_UVC(dx, dy, np.where(S.active.flatten(), 1.0, 0.3))
            else:
                q.set_UVC(dx, dy)
            ax.set_title(f't={P.initial+n*P.interval}')
            return q,
    else:
        fig, ax = plt.subplots(1,1)
        ax.set_xlim(0, P.shape[0])
        ax.set_ylim(0, P.shape[1])
        ax.set_title(f't={P.initial}')
        ax.set_aspect('equal')
        L = P.shape[0]*P.shape[1]
        q = ax.quiver( np.arange(L)//P.shape[1]+0.5, np.arange(L)%P.shape[1]+0.5, np.ones(L), np.zeros(L), pivot='mid',cmap=plt.get_cmap('Greys'),norm=Normalize(0,1) )
        # Update function
        def update(ns):
            n, S = ns
            dx,dy = P.dirs[S.directions.flatten()].T
            if args.active:
                q.set_UVC(dx, dy, np.where(S.active.flatten(), 1.0, 0.3))
            else:
                q.set_UVC(dx, dy)
            ax.set_title(f't={P.initial+n*P.interval}')
            return q,

    ani = FuncAnimation(fig, update, frames=enumerate(block), interval=50, blit=False, repeat=False, save_count=sys.maxsize)
    if args.fig is None:
        plt.show()
    else:
        # FIXME: at the moment, frames don't seem to get cleared; internet is silent on why
        ani.save(args.fig)

def blocks(args):
    """
    Movie of coloured blocks
    """

    block = process_block(args.file)
    # parameters will set the scene
    P = next(block)
    if(len(P.shape) > 2): raise RuntimeError("Can't (yet) create movies when d>2")

    if len(P.shape) == 1:
        # Determine aspect ratio
        aspect = max(1,P.shape[0]/20)
        fig, ax = plt.subplots(1,1, figsize=(6,1))
        ax.set_xlim(0, P.shape[0])
        ax.set_ylim(0, 1)
        ax.set_title(f't={P.initial}')
        im = ax.imshow(np.zeros((1,)+P.shape), extent=(0,P.shape[0],0,1), cmap=plt.get_cmap('Greys'), norm=Normalize(0,1), aspect=aspect)
        # Update function
        def update(ns):
            n, S = ns
            imdata = np.zeros(P.shape)
            imdata[np.where(S.occupied)] = 1.0
            if args.active:
                imdata[np.where(S.occupied ^ S.active)] = 0.2
            im.set_data(imdata.reshape((1,)+P.shape))
            ax.set_title(f't={P.initial+n*P.interval}')
            return im,
    else:
        fig, ax = plt.subplots(1,1)
        ax.set_xlim(0, P.shape[0])
        ax.set_ylim(0, P.shape[1])
        ax.set_aspect('equal')
        im = ax.imshow(np.zeros(P.shape), extent=(0,P.shape[0],0,P.shape[1]), cmap=plt.get_cmap('Greys'), norm=Normalize(0,1))
        # Update function
        def update(ns):
            n, S = ns
            imdata = np.zeros(P.shape)
            imdata[np.where(S.occupied)] = 1.0
            if args.active:
                imdata[np.where(S.occupied ^ S.active)] = 0.2
            im.set_data(imdata.T)
            ax.set_title(f't={P.initial+n*P.interval}')
            return im,

    pos = ax.get_position()

    ani = FuncAnimation(fig, update, frames=enumerate(block), interval=50, blit=False, repeat=False, save_count=sys.maxsize)

    if args.fig is None:
        plt.show()
    else:
        ani.save(args.fig)
