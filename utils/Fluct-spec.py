import numpy as np
import sys
import h5py
from scipy.interpolate import griddata
from scipy.fft import fft2, fftshift


def read_pos(filename):
    f = h5py.File(filename, 'r')
    pos = f["pos"][()]
    dim = int(len(pos) / 3)
    pos = pos.reshape(dim, 3)
    return pos


def radial_average(data, center=None):
    y, x = np.indices(data.shape)
    if center is None:
        center = np.array([data.shape[1] // 2, data.shape[0] // 2])
    r = np.hypot(x - center[0], y - center[1])
    r = r.astype(int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    return tbin / nr


snap_file = sys.argv[1]

N = 256
L = 2 * np.pi
dx = L / N
xx, yy = np.meshgrid(np.arange(0, L, dx), np.arange(0, L, dx))
fac = N**4 / (4 * np.pi**2)

pos = read_pos(snap_file)
F = griddata((pos[:, 0], pos[:, 1]), pos[:, 2], (xx, yy), method='linear')
F[np.isnan(F)] = 0
ft_data = fftshift(fft2(F))
power_spectrum = np.abs(ft_data)**2
spec = radial_average(power_spectrum) / fac
kk = np.linspace(1, len(spec), len(spec))
out_file = snap_file.replace("snap", "spec").replace("h5", "dat")
print(out_file)
np.savetxt(out_file, np.vstack([kk, spec]).T, fmt="%08e")
