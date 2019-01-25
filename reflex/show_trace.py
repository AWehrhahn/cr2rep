#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def compare(fname_trace, fname_img=None, fname_spec=None):
    """ compare img and trace """
    trace = fits.open(fname_trace)
    if fname_img:
        img = fits.open(fname_img)
        linecol = 'w'
    else: linecol= 'k'
    if fname_spec:
        spec = fits.open(fname_spec)

    X = np.arange(2048)
    FIG = plt.figure(figsize=(10, 3.5))

    for i in [1, 2, 3]:
        ax = FIG.add_subplot(1, 3, i)
        ax.set_xticks([])
        ax.set_yticks([])

        try: tdata = trace[i].data
        except:
            print('extension %s is missing, skipping.' % i)
            continue
        if tdata is None:
            print('Data for CHIP%s is empty, skipping.' % i)
            continue

        if fname_img:
            imgdata = img[i].data
            ax.imshow(imgdata, origin='lower',vmax=imgdata.max()*0.2)

        for t in tdata:
            upper = t['Upper']
            lower = t['Lower']
            alla = t['All']
            slitfrac = t['SlitFraction']
            order = t['Order']

            pol = np.polyval(t['Upper'][::-1], X)
            ax.plot(X, pol, ':'+linecol)

            pol = np.polyval(t['Lower'][::-1], X)
            ax.plot(X, pol, ':'+linecol)

            pol = np.polyval(t['All'][::-1], X)
            ax.plot(X, pol, '--'+linecol)
            if np.isnan(pol[1024]):
                continue
            #ax.text(1024, pol[1024], 'order: %s\ntrace: %s\nslitfrac: %.2f %.2f %.2f' % (t['order'], t['TraceNb'], *t['SlitFraction']), color=linecol, horizontalalignment='center',
            #        verticalalignment='center', size=11)


    FIG.tight_layout(pad=0.02)
    #plt.show()
    figname = fname_trace.replace('.fits','.png')
    plt.savefig(figname,dpi=120)


if __name__ == '__main__':
    fname_trace, fname_img, fname_spec = (None,)*3
    if len(sys.argv) > 1:
        fname_trace = sys.argv[1]
    if len(sys.argv) > 2:
        fname_img = sys.argv[2]
    if len(sys.argv) > 3:
        fname_spec = sys.argv[3]

    compare(fname_trace, fname_img, fname_spec)


