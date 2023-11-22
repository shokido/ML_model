import numpy as np
from ml_param import *
from ml_util import *
from nnf_sub import *
from solve_diag import *
import datetime as dt
import time
# Constants and parameters
fname_out = "out_case_nn.txt"
dt_start=dt.datetime(1900,1,1,0,0,0)
dt_end=dt.datetime(1900,1,11,0,0,0)
dz = 5.0
lat = 45.0
ntime = 6 * 24 * 20  # Adjust the number of time steps as needed
dt = 600.0; dt_output=60.0*60.0*24.0

bottom_depth = 1000.0
nu = 1.0e-6
nu_t = 1.0e-7
nu_s = 1.0e-7
beta1 = 0.6
beta2 = 20.0
r_long = 0.62

start = time.time()
ntime=int(((dt_end-dt_start).seconds+60*60*24*(dt_end-dt_start).days)/dt)
ntime_output=int(dt_output/dt)
coriolis = 2.0 * omega * np.sin(lat * np.pi / 180.0)
nz_rho=int(bottom_depth/dz)
nz_q=nz_rho+1
z_q=np.arange(-1*bottom_depth,dz,dz)
lev_q=-1*z_q
z_rho=np.arange(-1*bottom_depth+0.5*dz,0.5*dz,dz)
lev_rho=-1*z_rho

# Arrays initialization
temp_1d = np.zeros(nz_rho)
salt_1d = np.zeros(nz_rho)
u_1d = np.zeros(nz_rho)
v_1d = np.zeros(nz_rho)
temp_next_1d = np.zeros(nz_rho)
salt_next_1d = np.zeros(nz_rho)
u_next_1d = np.zeros(nz_rho)
v_next_1d = np.zeros(nz_rho)
temp_correct_1d = np.zeros(nz_rho)
salt_correct_1d = np.zeros(nz_rho)
u_correct_1d = np.zeros(nz_rho)
v_correct_1d = np.zeros(nz_rho)
bvf_1d = np.zeros(nz_rho - 1)
shear_1d = np.zeros(nz_rho - 1)
km_1d = np.zeros(nz_rho - 1)
kt_1d = np.zeros(nz_rho - 1)
ks_1d = np.zeros(nz_rho - 1)
penet_1d = np.zeros(nz_rho + 1)
qq_1d = np.zeros(nz_rho + 1)
l_1d = np.zeros(nz_rho + 1)
qq_next_1d = np.zeros(nz_rho + 1)
l_next_1d = np.zeros(nz_rho + 1)
kq_1d = np.zeros(nz_rho)


# Open output file for writing
with open(fname_out, "w") as output_file:
    output_file.write("z_rho\n")
    output_file.writelines([f"{z}\n" for z in z_rho])

# Set initial conditions for temperature, salinity, and velocity
for iz in range(nz_rho):
    if lev_rho[iz] <= 20.0:
        temp_1d[iz] = 15.0
    elif lev_rho[iz] <= 50.0:
        temp_1d[iz] = 15.0 + (9.0 - 15.0) * (lev_rho[iz] - 20.0) / (50.0 - 20.0)
    else:
        temp_1d[iz] = 9.0
    salt_1d[iz] = 35.0
    u_1d[iz] = 0.0
    v_1d[iz] = 0.0
# Initialize turbulence-related variables
qq_1d, l_1d = initialize_turb(nz_q, z_q)
# Main time-stepping loop
for itime in range(1, ntime + 1):
    hflx_nosolar = -20.0
    hflx_solar = 0.0
    sflx = 0.0
    uflx = 0.0
    vflx = 0.0
    #   Case 2
    #hflx_solar=0.0
    #sflx=0.0
    #hflx_nosolar=-70*np.sin(2.0*np.pi*(dt*itime/(60*60*24)/5))
    #uflx=0.2*np.sin(2.0*np.pi*(dt*itime/(60*60*24)/5))

    # Calculate shear and stratification
    bvf_1d = cal_bvf(z_rho, temp_1d, salt_1d)
    shear_1d = cal_shear(z_rho, u_1d, v_1d)
    penet_1d=cal_swheat( z_q, hflx_solar, beta1, beta2, r_long)
    km_1d, kt_1d, ks_1d=mynnf25_vdiff(nz_rho, bvf_1d, shear_1d, qq_1d, l_1d)

    km_1d=km_1d+nu
    kt_1d=kt_1d+nu_t
    ks_1d=ks_1d+nu_s
    kq_1d=cal_kq(nz_rho, km_1d)
    temp_next_1d=cal_next_theta(dt, nz_rho, z_rho, z_q, temp_1d, kt_1d, penet_1d, hflx_nosolar, temp_correct_1d)
    salt_next_1d=cal_next_salt(dt, nz_rho, z_rho, z_q, salt_1d, ks_1d, sflx, salt_correct_1d)
    u_next_1d,v_next_1d=cal_next_uv(dt, nz_rho, z_rho, z_q, u_1d, v_1d, km_1d, coriolis, uflx, vflx, u_correct_1d, v_correct_1d)
    qq_next_1d=cal_next_qq(dt, nz_rho, z_rho, z_q, qq_1d, bvf_1d, shear_1d, l_1d, kq_1d, km_1d, kt_1d, uflx, vflx, B_1)
#    l_next_1d=diag_l_mynnf(nz_rho, z_q, bvf_1d, qq_1d, uflx, vflx, hflx_solar, hflx_nosolar, sflx)
    l_next_1d=diag_l_mynnf(nz_rho, z_q, bvf_1d, qq_1d, uflx, vflx, hflx_solar, hflx_nosolar, sflx)
    
    # Update values
    temp_1d[:] = temp_next_1d[:]
    salt_1d[:] = salt_next_1d[:]
    u_1d[:] = u_next_1d[:]
    v_1d[:] = v_next_1d[:]
    qq_1d[:] = qq_next_1d[:]
    l_1d[:] = l_next_1d[:]
    
    
    if itime % (ntime_output) == 0:
        print(itime,temp_next_1d[nz_rho - 1])
        with open(fname_out, "a") as output_file:
            output_file.write(f"{itime * dt}\n")
            output_file.write("Temperature:\n")
            output_file.writelines([f"{temp}\n" for temp in temp_1d])
            output_file.write("Salinity:\n")
            output_file.writelines([f"{salt}\n" for salt in salt_1d])
            output_file.write("U Velocity:\n")
            output_file.writelines([f"{u}\n" for u in u_1d])
            output_file.write("V Velocity:\n")
            output_file.writelines([f"{v}\n" for v in v_1d])

print("Simulation completed.")
end = time.time()
print(end-start)
