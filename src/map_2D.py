import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import math
import sys

def read_map_data(filename):
    data = np.loadtxt(filename)

    xpos_ = data[:,0]
    ypos_ = data[:,1]
    cdm_dens_ = data[:,2]
    gas_dens_ = data[:,3]
    tmpr_ = data[:,4]
    fHI_ = data[:,5]

    nx=math.isqrt(data.shape[0])
    ny=math.isqrt(data.shape[0])    

    xpos = xpos_.reshape(nx,ny)
    ypos = ypos_.reshape(nx,ny)
    cdm_dens = cdm_dens_.reshape(nx,ny)
    gas_dens = gas_dens_.reshape(nx,ny)
    tmpr = tmpr_.reshape(nx,ny)
    fHI = fHI_.reshape(nx,ny)

    return xpos, ypos, cdm_dens, gas_dens, tmpr, fHI


def draw_map_2d(val, output_file, logscale):

    fig = plt.figure(dpi=300)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    ax = fig.add_axes((0.1, 0.1, 0.9, 0.9))
    ax.set_xlabel('x/L')
    ax.set_ylabel('y/L')

    if logscale:
        max_log10 = 10**math.floor(math.log10(val.max()))
        min_log10 = 10**math.ceil(math.log10(val.min()+1.0e-30))

#        max_log10 = min(max_log10, 10**2)
#        min_log10 = max(min_log10, 10**-1)

        log_norm = LogNorm(vmin=min_log10, vmax=max_log10)

        im=ax.imshow(val+1.0e-10,origin='lower',norm=log_norm,extent=(0.0, 1.0, 0.0, 1.0),interpolation='none')
    else:
        im=ax.imshow(val+1.0e-10,origin='lower',extent=(0.0, 1.0, 0.0, 1.0),interpolation='none')

    #ax.quiver(velx, vely)

    fig.colorbar(im, ax=ax)
    plt.savefig(output_file, bbox_inches='tight')

xx, yy, cdm_dens, gas_dens, tmpr, fHI = read_map_data(sys.argv[1])

print(gas_dens.min(), gas_dens.max())
draw_map_2d(gas_dens, 'gas_dens_map.pdf', True)
draw_map_2d(tmpr, 'gas_tmpr_map.pdf', True)
draw_map_2d(1.0-fHI, 'gas_fHI_map.pdf', False)



