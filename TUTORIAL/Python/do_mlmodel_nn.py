import numpy as np
from ml_param import *
from ml_util import *
from nnf_sub import *
from solve_diag import *
import datetime as dtime
import time
# Constants and parameters
fname_out = "out_case_nn.txt"
dt_start=dtime.datetime(1900,1,1,0,0,0)
dt_end=dtime.datetime(1900,1,21,0,0,0)
dz = 5.0
lat = 45.0
ntime = 6 * 24 * 20  # Adjust the number of time steps as needed
dt = 240.0; dt_output=60.0*60.0*24.0

bottom_depth = 1000.0
NB=-1 # Reaching bottom
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
lev_q=np.arange(0,bottom_depth+dz,dz)
z_q=-1*lev_q
z_rho=0.5*(z_q[0:nz_q-1]+z_q[1:nz_q])
lev_rho=-1*z_rho

# Arrays initialization
temp_1d = np.zeros(nz_rho);salt_1d = np.zeros(nz_rho)
u_1d = np.zeros(nz_rho);v_1d = np.zeros(nz_rho)
temp_next_1d = np.zeros(nz_rho);salt_next_1d = np.zeros(nz_rho)
u_next_1d = np.zeros(nz_rho);v_next_1d = np.zeros(nz_rho)
temp_past_1d = np.zeros(nz_rho);salt_past_1d = np.zeros(nz_rho)
u_past_1d = np.zeros(nz_rho);v_past_1d = np.zeros(nz_rho)

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
qq_1d = np.zeros(nz_rho + 1);l_1d = np.zeros(nz_rho + 1)
qq_next_1d = np.zeros(nz_rho + 1);l_next_1d = np.zeros(nz_rho + 1)
qq_past_1d = np.zeros(nz_rho + 1);l_past_1d = np.zeros(nz_rho + 1)
kq_1d = np.zeros(nz_rho)

nstep_out=0
for itime in range(1, ntime + 1):
    if itime % (ntime_output) == 0:
        nstep_out+=1
# Open output file for writing
with open(fname_out, "w") as output_file:
    output_file.write(str(nstep_out)+"\n")
    out_list=[f"{ivar} " for ivar in lev_rho];out_list.append("\n")
    output_file.writelines(out_list)

# Set initial conditions for temperature, salinity, and velocity
for iz in range(nz_rho):
    if lev_rho[iz] <= 20.0:
        temp_1d[iz] = 15.0
    elif lev_rho[iz] <= 50.0:
        temp_1d[iz] = 15.0 + (12.0 - 15.0) * (lev_rho[iz] - 20.0) / (50.0 - 20.0)
    else:
        temp_1d[iz] = 12.0
    salt_1d[iz] = 35.0
    u_1d[iz] = 0.0
    v_1d[iz] = 0.0
# Initialize turbulence-related variables
qq_1d, l_1d = initialize_turb(nz_q, z_q)
temp_past_1d=np.copy(temp_1d);salt_past_1d=np.copy(salt_1d)
u_past_1d=np.copy(u_1d);v_past_1d=np.copy(v_1d)
qq_past_1d=np.copy(qq_1d);qq_past_1d=np.copy(qq_1d)
# Main time-stepping loop
for itime in range(1, ntime + 1):
    dt_now=dt_start+dtime.timedelta(seconds=int(itime*dt))
    hflx_nosolar = 0.0;hflx_solar = 0.0
    sflx = 0.0
    uflx = 0.2;vflx = 0.0
    # #   Case 2
    #    hflx_solar=0.0;hflx_nosolar=-30.0
    #    sflx = 0.0
    #    uflx = 0.0;vflx = 0.0
    # #   Case 3
    #    hflx_solar=0.0;hflx_nosolar=30.0
    #    sflx = 0.0
    #    uflx = 0.0;vflx = 0.0
    # #   Case 4
    #    hflx_solar=0.0;hflx_nosolar=-30.0
    #    sflx = 0.0
    #    uflx = 0.2;vflx = 0.0
    # #   Case 5
    #    hflx_solar=0.0;hflx_nosolar=0.0
    #    uflx = 0.2;vflx = 0.0
    #  sflx=-2.5e-5 ! (E-P)
    # Calculate shear and stratification
    bvf_1d = cal_bvf(z_rho, temp_1d, salt_1d)
    shear_1d = cal_shear(z_rho, u_1d, v_1d)
    penet_1d=cal_swheat( z_q, hflx_solar, beta1, beta2, r_long)
    km_1d, kt_1d, ks_1d=mynnf25_vdiff(nz_rho, bvf_1d, shear_1d, qq_1d, l_1d)
    km_1d=km_1d+nu
    kt_1d=kt_1d+nu_t
    ks_1d=ks_1d+nu_s
    kq_1d=cal_kq(nz_rho, km_1d,NB=NB)
    # Time stepping
    temp_next_1d=cal_next_theta(2*dt, nz_rho, z_rho, z_q, temp_past_1d, kt_1d, penet_1d, hflx_nosolar, temp_correct_1d,NB=NB)
    salt_next_1d=cal_next_salt(2*dt, nz_rho, z_rho, z_q, salt_past_1d, ks_1d, sflx, salt_correct_1d,NB=NB)
    u_next_1d,v_next_1d=cal_next_uv(2*dt, nz_rho, z_rho, z_q, u_past_1d, v_past_1d, km_1d, coriolis, uflx, vflx, u_correct_1d, v_correct_1d,NB=NB)
    qq_next_1d=cal_next_qq(2*dt, nz_rho, z_rho, z_q, qq_past_1d, bvf_1d, shear_1d, l_1d, kq_1d, km_1d, kt_1d, uflx, vflx, B_1,NB=NB)
    l_next_1d=diag_l_mynnf(nz_rho, z_q, bvf_1d, qq_1d, uflx, vflx, hflx_solar, hflx_nosolar, sflx,NB=NB)
    # Update values
    temp_past_1d=np.copy(temp_1d);salt_past_1d=np.copy(salt_1d);
    u_past_1d=np.copy(u_1d);v_past_1d=np.copy(v_1d);qq_past_1d=np.copy(qq_1d)
    temp_1d=np.copy(temp_next_1d);salt_1d=np.copy(salt_next_1d);
    u_1d=np.copy(u_next_1d);v_1d=np.copy(v_next_1d);qq_1d=np.copy(qq_next_1d)
    l_1d[:] = np.copy(l_next_1d[:])

    if itime % (ntime_output) == 0:
        print(dt_now.year,dt_now.month,dt_now.day,temp_next_1d[0])
        with open(fname_out, "a") as output_file:
            output_file.write(str(dt_now.year)+" "+str(dt_now.month)+" "+str(dt_now.day)+" "+str(dt_now.hour)+" "+str(dt_now.minute)+" "+str(dt_now.second)+"\n")
            out_list=[f"{ivar} " for ivar in temp_1d];out_list.append("\n")
            output_file.writelines(out_list)
            out_list=[f"{ivar} " for ivar in salt_1d];out_list.append("\n")
            output_file.writelines(out_list)
            out_list=[f"{ivar} " for ivar in u_1d];out_list.append("\n")
            output_file.writelines(out_list)
            out_list=[f"{ivar} " for ivar in v_1d];out_list.append("\n")
            output_file.writelines(out_list)

print("Simulation completed.")
end = time.time()
print(end-start)
