import numpy as np
import matplotlib.pyplot as plt

def build_model(nx, nz, property, z):
    
    model = np.zeros((nz, nx), dtype = float)

    model[:int(z[0]),:] = np.ones((int(z[0]), nx)) * property[0]

    for depth in range(1, len(z)):
        
        layer = slice(int(z[depth - 1]), int(z[depth]))

        model[layer,:] = np.ones((int(z[depth] - z[depth-1]), nx)) * property[depth]

    return model

nx = 2001
nz = 201

dh = 10.0

z = np.array([500, 800, 1100, 1400, 1700, nz*dh]) // dh
vp = np.array([1500, 1650, 1800, 2000, 2200, 2500])
vs = vp / 1.7; vs[0] = 0.0
rho = 310 * vp ** 0.25; rho[0] = 1000

vp_model = build_model(nx, nz, vp, z)
vs_model = build_model(nx, nz, vs, z)
rho_model = build_model(nx, nz, rho, z)

plt.figure(1, figsize = (15, 8))

plt.subplot(311)
plt.imshow(vp_model)

plt.subplot(312)
plt.imshow(vs_model)

plt.subplot(313)
plt.imshow(rho_model)

plt.show()

folder = "../../../inputs/model/"

vp_model.flatten("F").astype("float32", order = "F").tofile(folder + f"layer_cake_vp_{nz}x{nx}_{dh:.0f}m.bin")
vs_model.flatten("F").astype("float32", order = "F").tofile(folder + f"layer_cake_vs_{nz}x{nx}_{dh:.0f}m.bin")
rho_model.flatten("F").astype("float32", order = "F").tofile(folder + f"layer_cake_rho_{nz}x{nx}_{dh:.0f}m.bin")
