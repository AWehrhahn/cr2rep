#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def compare(fname_trace, fname_img=None, fname_spec=None):
    """ compare img and trace """
    trace = fits.open(fname_trace)
    if fname_img:
        img = fits.open(fname_img)
        linecol = "w"
    else:
        linecol = "k"
    if fname_spec:
        spec = fits.open(fname_spec)

    X = np.arange(2048)
    FIG = plt.figure(figsize=(10, 3.5))

    for i in [1, 2, 3]:
        ax = FIG.add_subplot(1, 3, i)
        ax.set_xticks([])
        ax.set_yticks([])

        try:
            tdata = trace[i].data
        except:
            print("extension %s is missing, skipping." % i)
            continue
        if tdata is None:
            print("Data for CHIP%s is empty, skipping." % i)
            continue

        if fname_img:
            imgdata = img[i].data
            vmin, vmax = np.percentile(imgdata, (5, 95))
            ax.imshow(imgdata, origin="lower", vmin = vmin, vmax=vmax)

        for t in tdata:
            upper = t["Upper"]
            lower = t["Lower"]
            alla = t["All"]
            slitfrac = t["SlitFraction"]
            order = t["Order"]
            wave = t["Wavelength"]
            ca = t["SlitPolyA"]
            cb = t["SlitPolyB"]
            cc = t["SlitPolyC"]

            upper = np.polyval(upper[::-1], X)
            ax.plot(X, upper, ":" + linecol)

            lower = np.polyval(lower[::-1], X)
            ax.plot(X, lower, ":" + linecol)

            middle = np.polyval(alla[::-1], X)
            ax.plot(X, middle, "--" + linecol)

            i1 = tdata[tdata["order"] == order]["Slitfraction"][:, 1]
            i2 = tdata[tdata["order"] == order]["All"]
            coeff = [np.interp(0.5, i1, i2[:, k]) for k in range(i2.shape[1])]

            a, b, c = wave
            xw = np.polyval(wave[::-1], np.arange(2048))
            x = [(np.sqrt(-4*a*c + b**2 + 4 * c * w) - b) / (2*c) for w in np.arange(np.floor(xw.min()), np.ceil(xw.max()))]
            x = np.clip(np.asarray(x), 0, 2048)
            idx = x.astype(int)
            plt.vlines(x, lower[idx], upper[idx], colors="r")


            for i in range(10, 2048, 100):
                ew = [int(middle[i] - lower[i]), int(upper[i] - middle[i])]
                x = np.zeros(ew[0] + ew[1] + 1)
                y = np.arange(-ew[0], ew[1] + 1).astype(float)

                # Evaluate the curvature polynomials to get coefficients
                a = np.polyval(ca[::-1], i)
                b = np.polyval(cb[::-1], i)
                c = np.polyval(cc[::-1], i)
                yc = np.polyval(coeff[::-1], i)

                # Shift polynomials to the local frame of reference
                a = a - i + yc * b + yc * yc * c
                b += 2 * yc * c

                for j, yt in enumerate(y):
                    x[j] = i + yt * b + yt ** 2 * c
                y += middle[i]
                plt.plot(x, y, "--r")

            if np.isnan(middle[1024]):
                continue
            ax.text(
                1024,
                middle[1024],
                "order: %s\ntrace: %s" % (t["order"], t["TraceNb"]),
                color=linecol,
                horizontalalignment="center",
                verticalalignment="center",
                size=7,
            )

    FIG.tight_layout(pad=0.02)
    plt.show()
    figname = fname_trace.replace(".fits", ".png")
    plt.savefig(figname, dpi=120)


if __name__ == "__main__":
    compare(*sys.argv[1:])