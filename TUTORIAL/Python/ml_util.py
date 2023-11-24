import numpy as np
from ml_param import *
# Constants and parameters
kappa = 0.41
q2_init = 1.0e-8
tiny = 1.0e-12

def cal_pdens(s, t):
    # Calculate seawater potential density (rho) as a function of salinity (s) and potential temperature (t)
    rho_w = (999.842594 + t * 6.793952e-2 - t**2 * 9.095290e-3 +
             1.001685e-04 * t**3 - 1.120083e-6 * t**4 + 6.536332e-9 * t**5)
    
    ret = rho_w + (0.824493 - 4.0899e-3 * t + 7.6438e-5 * t**2 - 8.2467e-7 * t**3 +
                   5.3875e-9 * t**4) * s + (-5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t**2) * s**1.5 + 4.8314e-4 * s**2
    return ret

def cal_alpha(s, t):
    # Calculate the thermal expansion coefficient (alpha) of seawater as a function of salinity (s) and temperature (t)
    dt = 0.01
    rho_1 = cal_pdens(s, t)
    rho_2 = cal_pdens(s, t + dt)
    alpha = (rho_2 - rho_1) / (rho_1 * dt)
    return alpha

def cal_beta(s, t):
    # Calculate the salinity contraction coefficient (beta) of seawater as a function of salinity (s) and temperature (t)
    ds = 0.01
    rho_1 = cal_pdens(s + ds, t)
    rho_2 = cal_pdens(s, t)
    beta = (rho_1 - rho_2) / (rho_1 * ds)
    return beta

def cal_bvf(z_rho, temp, salt,NB=-1):
    # Calculate the buoyancy frequency (N2) of seawater as a function of vertical grid information (z_rho),
    # temperature (temp), and salinity (salt)
    nz = len(z_rho)
    if (NB==-1):
        NB=nz
#    diff_pdens=(cal_pdens(salt[0:(NB-1)],temp[0:(NB-1)])-cal_pdens(salt[1:NB],temp[1:NB]))/cal_pdens(salt[0:(NB-1)],temp[0:(NB-1)])
    diff_pdens=cal_pdens(salt[0:(NB-1)],temp[0:(NB-1)])-cal_pdens(salt[1:NB],temp[1:NB])
    sum_pdens=0.5*(cal_pdens(salt[0:(NB-1)],temp[0:(NB-1)])+cal_pdens(salt[1:NB],temp[1:NB]))
    bv=-1*g*diff_pdens/(sum_pdens*(z_rho[0:(NB-1)]-z_rho[1:NB]))
    return bv

def cal_shear(z_rho, u, v,NB=-1):
    # Calculate the shear of ocean currents as a function of vertical grid information (z_rho),
    # and horizontal velocity components (u and v)
    nz = len(z_rho)
    if (NB==-1):
        NB=nz
    shear = np.zeros(NB - 1)
    diff_z=z_rho[0:(NB-1)]-z_rho[1:NB]
    uz=u[0:(NB-1)]-u[1:NB]
    vz=v[0:(NB-1)]-v[1:NB]
    shear=((uz/diff_z)**2 +(vz/diff_z)**2+tiny)    
    return shear

def initialize_turb(N, z):
    # Initialize turbulence-related variables (q2 and l) based on vertical grid information (z)
    H = 20.0
    q2 = q2_init * np.exp(z / H)
    l = kappa * np.abs(z)
    return q2, l

def cal_penet(level, sw, beta1, beta2, r_long):
    # Calculate the penetration of solar radiation into the water column
    z1 = abs(level)
    z1b1 = min(z1 / beta1, 100.0)
    z1b2 = min(z1 / beta2, 100.0)
    
    ret = sw * (r_long * np.exp(-z1b1) + (1.0 - r_long) * np.exp(-z1b2))
    return ret

def cal_swheat(level, heat_sl, beta1, beta2, r_long):
    # Calculate heat absorption in the water column based on the penetration of solar radiation
    nz = len(level)
    heat = np.zeros(nz)
    heat=[cal_penet(i,heat_sl, beta1, beta2, r_long) for i in level]
    heat=np.asarray(heat)
    return heat

def cal_mo_inv(sst, sss, tau_x, tau_y, sw, nsw, salflux):
    # Calculate the Monin-Obukhov stability parameter (Mo_inv)
    alpha = cal_alpha(sss, sst)
    beta = cal_beta(sss, sst)
    bf = -g * (alpha * (sw + nsw) / (rho * cpw) + beta * salflux)
    
    tau_x = np.abs(tau_x)
    tau_y = np.abs(tau_y)
    u_star2 = np.sqrt(tau_x**2 + tau_y**2) / rho
    u_star3 = np.sqrt(u_star2)**3
    
    u_star3 = max(u_star3, 1.0e-10)
    
    mo_inv = kappa * bf / u_star3
    return mo_inv
