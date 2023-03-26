import sys 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

def catch_parameter(target, filename):
    file = open(filename,'r')
    for line in file.readlines():
        if line[0] != '#':
            splitted = line.split()
            if len(splitted) != 0:    
                if splitted[0] == target:
                    if splitted[2] == '[':
                        output = np.array([])
                        for s in splitted[3:]:
                            if s not in ']': 
                                if s[-1] == ',':
                                    output = np.append(output, float(s[:-1]))        
                                else: 
                                    output = np.append(output, float(s))        
                            else:
                                break

                        return output
                    else:
                        return splitted[2]

def build_model(nx, nz, property, z):
    model = np.zeros((nz, nx), dtype = float)
    model[:int(z[0]),:] = np.ones((int(z[0]), nx)) * property[0]

    for depth in range(1, len(z)):
        layer = slice(int(z[depth - 1]), int(z[depth]))
        model[layer,:] = np.ones((int(z[depth] - z[depth-1]), nx)) * property[depth]

    return model

file = sys.argv[1]

nx = int(catch_parameter("nx", file))
nz = int(catch_parameter("nz", file))

dx = float(catch_parameter("dx", file))
dz = float(catch_parameter("dz", file))

z = catch_parameter("z", file)
vp = catch_parameter("vp", file)
vs = catch_parameter("vs", file)
rho = catch_parameter("rho", file)

z = np.append(z, nz*dz)

z /= dz

vp_model = build_model(nx, nz, vp, z)
vs_model = build_model(nx, nz, vs, z)
rho_model = build_model(nx, nz, rho, z)

prop_min = np.min(np.append(np.append(vp, vs), rho)) 
prop_max = np.max(np.append(np.append(vp, vs), rho))  

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc * dx, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc * dz, dtype = int)

ncbar = 5

labels = ["P velocity [m/s]", "S velocity [m/s]", "Density [kg/m³]"]

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (15,8))

ax[0].imshow(vp_model, aspect = "auto", cmap = "Greys", vmin = prop_min, vmax = prop_max)
ax[1].imshow(vs_model, aspect = "auto", cmap = "Greys", vmin = prop_min, vmax = prop_max)
ax[2].imshow(rho_model, aspect = "auto", cmap = "Greys", vmin = prop_min, vmax = prop_max)

for k in range(len(ax)):
    ax[k].set_xlabel("Distance [m]", fontsize = 15)
    ax[k].set_ylabel("Depth [m]", fontsize = 15) 

    ax[k].set_xticks(xloc)
    ax[k].set_xticklabels(xlab)

    ax[k].set_yticks(zloc)
    ax[k].set_yticklabels(zlab)

    cmap = mpl.colormaps["Greys"]
    norm = mpl.colors.Normalize(prop_min, prop_max)
    divider = make_axes_locatable(ax[k])
    cax = divider.append_axes("right", size="1%", pad=0.5)

    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, ticks = np.linspace(prop_min, prop_max, ncbar), orientation = "vertical")
    cbar.set_label(labels[k], fontsize = 12)

plt.tight_layout()
plt.show()

vp_file = catch_parameter("vp_file", file)
vs_file = catch_parameter("vs_file", file)
rho_file = catch_parameter("rho_file", file)

vp_model.flatten("F").astype("float32", order = "F").tofile(vp_file)
vs_model.flatten("F").astype("float32", order = "F").tofile(vs_file)
rho_model.flatten("F").astype("float32", order = "F").tofile(rho_file)
