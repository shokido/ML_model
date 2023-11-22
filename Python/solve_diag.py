import numpy as np
from ml_param import *
from scipy.linalg import solve_banded
# Constants
qq_min = 1.0e-8

# Subroutine for solving a tridiagonal system
def solve_tri_implicit(A_in, B_in, C_in, D_in):
    N=np.shape(A_in)[0]
    ab = np.zeros((3, N), dtype=complex)
    ab[0, 1:] = A_in[1:N]
    ab[1, :] = B_in[0:N]
    ab[2, :-1] = C_in[0:(N-1)]
    X = solve_banded((1, 1), ab, D_in)
    return X
# Subroutine for calculating the next temperature
def cal_next_theta(dt, N, z_rho, z_q, theta, kt, heat_penet, heat_no_solar, temp_correct):
    # z_q[i]---------------z_q[i+1]
    #          z_rho[i]---------------z_rho[i+1]
    # kt[i-1]-------kt[i]
    A = np.zeros(N);B = np.zeros(N);C = np.zeros(N);D = np.zeros(N)
    A[1:N]=-1*dt*np.copy(kt[0:(N-1)])/((z_q[2:(N+1)]-z_q[1:N])*(z_rho[1:N]-z_rho[0:(N-1)]))
    C[0:N-1]=-1*dt*np.copy(kt[0:(N-1)])/((z_q[1:N]-z_q[0:(N-1)])*(z_rho[1:N]-z_rho[0:(N-1)]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]
    D[1:(N-1)]=theta[1:(N-1)]+dt*(heat_penet[1:(N-1)]-heat_penet[0:(N-2)])/(rho*cpw*(z_q[2:N]-z_q[1:(N-1)]))+dt*temp_correct[1:(N-1)]
    # Bottom boundary conditions
    A[0]=0.0
    B[0]=1.0-C[0]
    D[0]=theta[0]
    # Top boundary conditions
    B[N-1]=1.0-A[N-1]
    C[N-1]=0
    D[N-1]=theta[N-1]+dt*(heat_penet[N-1]-heat_penet[N-2]+heat_no_solar)/(rho*cpw*(z_q[N]-z_q[N-1]))+dt*temp_correct[N-1]
    # Solve tridiagonal equations
#    theta_next=solve_tri_implicit(A, B, C, D)
    theta_next=solve_tri_implicit(A, B, C, D)
    return(np.real(theta_next))
# Subroutine for calculating the next salinity
def cal_next_salt(dt, N, z_rho, z_q, salt, ks, ssflux, salt_correct):
    salt_next=np.zeros(N)
    A = np.zeros(N);B = np.zeros(N);C = np.zeros(N);D = np.zeros(N)
    A[1:N]=-1*dt*np.copy(ks[0:(N-1)])/((z_q[2:(N+1)]-z_q[1:N])*(z_rho[1:N]-z_rho[0:(N-1)]))
    C[0:N-1]=-1*dt*np.copy(ks[0:(N-1)])/((z_q[1:N]-z_q[0:(N-1)])*(z_rho[1:N]-z_rho[0:(N-1)]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]
    D[1:(N-1)]=salt[1:(N-1)]+dt*salt_correct[1:(N-1)]
    # Bottom boundary conditions
    A[0]=0.0
    B[0]=1.0-C[0]
    D[0]=salt[0]
    # Top boundary conditions
    B[N-1]=1.0-A[N-1]
    C[N-1]=0
    D[N-1]=salt[N-1]+dt*ssflux/(z_q[N]-z_q[N-1])+dt*salt_correct[N-1]
#    salt_next=solve_tri_implicit(A, B, C, D)
    salt_next=solve_tri_implicit(A, B, C, D)
    return(np.real(salt_next))
# Subroutine for calculating the next velocity components (u and v)
def cal_next_uv(dt, N, z_rho, z_q, u, v, km, f, tau_x, tau_y, u_correct, v_correct):
    u_next=np.zeros(N);v_next=np.zeros(N)
    A = np.zeros(N,dtype=complex);B = np.zeros(N,dtype=complex)
    C = np.zeros(N,dtype=complex);D = np.zeros(N,dtype=complex)
    uv=u + 1j * v
    uv_correct=u_correct + 1j * v_correct
    A[1:N]=-1*dt*np.copy(km[0:(N-1)])/((z_q[2:(N+1)]-z_q[1:N])*(z_rho[1:N]-z_rho[0:(N-1)]))
    C[0:N-1]=-1*dt*np.copy(km[0:(N-1)])/((z_q[1:N]-z_q[0:(N-1)])*(z_rho[1:N]-z_rho[0:(N-1)]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]+np.ones(N-2)*f*dt*1j
    D[1:(N-1)]=uv[1:(N-1)]+dt*uv_correct[1:(N-1)]
    # Bottom boundary conditions
    A[0]=0.0
    B[0]=1.0-C[0]
    D[0]=uv[0]
    # Top boundary conditions
    B[N-1]=1.0-A[N-1]+dt*f*1j
    C[N-1]=0
    D[N-1]=uv[N-1]+dt*(tau_x+1j*tau_y)/(rho*(z_q[N]-z_q[N-1]))+dt*uv_correct[N-1]
#    uv_next=solve_tri_implicit(A, B, C, D)
    uv_next=solve_tri_implicit(A, B, C, D)
    u_next = np.real(uv_next)
    v_next = np.imag(uv_next)
    return(u_next,v_next)
# # Subroutine for calculating the next turbulent kinetic energy (TKE or qq)
def cal_next_qq(dt, N, z_rho, z_q, qq, bvf, shear2, l, kq, km, kt, tau_x, tau_y, B_1):
    qq_next=np.zeros(N+1)
    A = np.zeros(N+1);B = np.zeros(N+1);C = np.zeros(N+1);D = np.zeros(N+1)
    SP = np.zeros(N+1);BP = np.zeros(N+1);DS = np.zeros(N+1)
    
    A[1:N]=-1*dt*np.copy(kq[0:(N-1)])/((z_rho[1:N]-z_rho[0:(N-1)])*(z_q[1:N]-z_q[0:(N-1)]))
    C[1:N]=-1*dt*np.copy(kq[1:N])/((z_rho[1:N]-z_rho[0:(N-1)])*(z_q[2:(N+1)]-z_q[1:N]))    
    B[1:N]=1.0-A[1:N]-C[1:N]
    SP[1:N]=km[0:(N-1)]*shear2[0:(N-1)]
    BP[1:N]=-1*kt[0:(N-1)]*bvf[0:(N-1)]
    DS[1:N]=(qq[1:(N)]**1.5)/(B_1*l[1:N])
    D[1:N]=qq[1:N]+2*dt*(SP[1:N]+BP[1:N]-DS[1:N])

    A[0]=0.0;B[0]=1.0;C[0]=0;D[0]=0
    # Top boundary conditions
    tau_ts = np.sqrt(tau_x**2 + tau_y**2) / rho
    A[N]=0
    B[N]=1.0
    C[N] = 0.0
    D[N] = (B_1**(2.0 / 3.0)) * tau_ts
    qq_next=solve_tri_implicit(A, B, C, D)
    qq_next = np.maximum(qq_next, qq_min)
    return(np.real(qq_next))
