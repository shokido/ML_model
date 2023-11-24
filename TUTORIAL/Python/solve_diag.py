import numpy as np
from ml_param import *
from scipy.linalg import solve_banded
# Constants
qq_min = 1.0e-8

# Subroutine for solving a tridiagonal system
def solve_tri_implicit(A_in, B_in, C_in, D_in):
    N=np.shape(A_in)[0]
    ab = np.zeros((3, N), dtype=complex)
    ab[0, 1:] = C_in[0:(N-1)]
    ab[1, :] = B_in[0:N]
    ab[2, :-1] = A_in[1:N]
    X = solve_banded((1, 1), ab, D_in)
    return X
# Subroutine for calculating the next temperature
def cal_next_theta(dt, N,z_rho, z_q, theta, kt, heat_penet, heat_no_solar, temp_correct,NB=-1):
    if (NB==-1):
        NB=N
    # z_q[i]---------------z_q[i+1]
    #          z_rho[i]---------------z_rho[i+1]
    # kt[i-1]--------------kt[i]
    A = np.zeros(N);B = np.zeros(N);C = np.zeros(N);D = np.zeros(N)
    theta_next = np.zeros(N)
    A[1:N]=-1*dt*np.copy(kt[0:(N-1)])/((z_q[1:N]-z_q[2:(N+1)])*(z_rho[0:(N-1)]-z_rho[1:N])) # Upper
    C[0:N-1]=-1*dt*np.copy(kt[0:(N-1)])/((z_q[0:(N-1)]-z_q[1:N])*(z_rho[0:(N-1)]-z_rho[1:N]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]
    D[1:(N-1)]=theta[1:(N-1)]+dt*(heat_penet[0:(N-2)]-heat_penet[1:(N-1)])/(rho*cpw*(z_q[1:(N-1)]-z_q[2:N]))+dt*temp_correct[1:(N-1)]
    # Top boundary conditions
    B[0]=1.0-C[0]
    A[0]=0
    D[0]=theta[0]+dt*(heat_penet[0]-heat_penet[1]+heat_no_solar)/(rho*cpw*(z_q[0]-z_q[1]))+dt*temp_correct[0]
    # Bottom boundary conditions
    B[NB-1]=1.0-A[NB-1]
    C[NB-1]=0.0
    D[NB-1]=theta[NB-1]+dt*temp_correct[NB-1]
    # Solve tridiagonal equations
#    theta_next=solve_tri_implicit(A, B, C, D)
    tmp=solve_tri_implicit(A[0:NB], B[0:NB], C[0:NB], D[0:NB])
    theta_next[0:NB]=np.real(tmp)
    if (NB < N):
        theta_next[NB:]=np.nan
    return(theta_next)
# Subroutine for calculating the next salinity
def cal_next_salt(dt, N, z_rho, z_q, salt, ks, ssflux, salt_correct,NB=-1):
    if (NB==-1):
        NB=N
    salt_next=np.zeros(N)
    A = np.zeros(N);B = np.zeros(N);C = np.zeros(N);D = np.zeros(N)
    salt_next = np.zeros(N)
    A[1:N]=-1*dt*np.copy(ks[0:(N-1)])/((z_q[1:N]-z_q[2:(N+1)])*(z_rho[0:(N-1)]-z_rho[1:N]))
    C[0:N-1]=-1*dt*np.copy(ks[0:(N-1)])/((z_q[0:(N-1)]-z_q[1:N])*(z_rho[0:(N-1)]-z_rho[1:N]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]
    D[1:(N-1)]=salt[1:(N-1)]+dt*salt_correct[1:(N-1)]
    # Top boundary conditions
    B[0]=1.0-C[0]
    A[0]=0
    D[0]=salt[0]+dt*ssflux/(z_q[0]-z_q[1])+dt*salt_correct[0]
    # Bottom boundary conditions
    B[NB-1]=1.0-A[NB-1]
    C[NB-1]=0.0
    D[NB-1]=salt[NB-1]
#    salt_next=solve_tri_implicit(A, B, C, D)
    tmp=solve_tri_implicit(A[0:NB], B[0:NB], C[0:NB], D[0:NB])
    salt_next[0:NB]=np.real(tmp)
    if (NB < N):
        salt_next[NB:]=np.nan
    return(salt_next)
# Subroutine for calculating the next velocity components (u and v)
def cal_next_uv(dt, N, z_rho, z_q, u, v, km, f, tau_x, tau_y, u_correct, v_correct,NB=-1):
    if (NB==-1):
        NB=N
    u_next=np.zeros(N);v_next=np.zeros(N)
    A = np.zeros(N,dtype=complex);B = np.zeros(N,dtype=complex)
    C = np.zeros(N,dtype=complex);D = np.zeros(N,dtype=complex)
    u_next = np.zeros(N);v_next = np.zeros(N)

    uv=u + 1j * v
    uv_correct=u_correct + 1j * v_correct
    A[1:N]=-1*dt*np.copy(km[0:(N-1)])/((z_q[1:N]-z_q[2:(N+1)])*(z_rho[0:(N-1)]-z_rho[1:N]))
    C[0:N-1]=-1*dt*np.copy(km[0:(N-1)])/((z_q[0:(N-1)]-z_q[1:N])*(z_rho[0:(N-1)]-z_rho[1:N]))
    B[1:(N-1)]=1.0-A[1:(N-1)]-C[1:(N-1)]+np.ones(N-2)*f*dt*1j
    D[1:(N-1)]=uv[1:(N-1)]+dt*uv_correct[1:(N-1)]
    # Top boundary conditions
    B[0]=1.0-C[0]+dt*f*1j
    A[0]=0
    D[0]=uv[0]+dt*(tau_x+1j*tau_y)/(rho*(z_q[0]-z_q[1]))+dt*uv_correct[0]
    # Bottom boundary conditions
#    B[NB-1]=1.0-A[NB-1]
#    C[NB-1]=0.0
#    D[NB-1]=uv[NB-1]
    B[NB-1]=1.0-A[NB-1]+dt*f*1j
    C[NB-1]=0.0
    D[NB-1]=uv[NB-1]+dt*uv_correct[NB-1]

#    uv_next=solve_tri_implicit(A, B, C, D)
    tmp=solve_tri_implicit(A[0:NB], B[0:NB], C[0:NB], D[0:NB])
    u_next[0:NB]=np.real(tmp)
    v_next[0:NB]=np.imag(tmp)
    if (NB < N):
        u_next[NB:]=np.nan
        v_next[NB:]=np.nan
    return(u_next,v_next)



# # Subroutine for calculating the next turbulent kinetic energy (TKE or qq)
def cal_next_qq(dt, N, z_rho, z_q, qq, bvf, shear2, l, kq, km, kt, tau_x, tau_y, B_1,NB=-1):
    if (NB==-1):
        NB=N
    qq_next=np.zeros(N+1)
    A = np.zeros(NB+1);B = np.zeros(NB+1);C = np.zeros(NB+1);D = np.zeros(NB+1)
    SP = np.zeros(NB+1);BP = np.zeros(NB+1);DS = np.zeros(NB+1)
    
    A[1:NB]=-1*dt*np.copy(kq[0:(NB-1)])/((z_rho[0:(NB-1)]-z_rho[1:NB])*(z_q[0:(NB-1)]-z_q[1:NB]))
    C[1:NB]=-1*dt*np.copy(kq[1:NB])/((z_rho[0:(NB-1)]-z_rho[1:NB])*(z_q[1:NB]-z_q[2:(NB+1)]))    
    B[1:NB]=1.0-A[1:NB]-C[1:NB]
    SP[1:NB]=km[0:(NB-1)]*shear2[0:(NB-1)]
    BP[1:NB]=-1*kt[0:(NB-1)]*bvf[0:(NB-1)]
    DS[1:NB]=(qq[1:NB]**1.5)/(B_1*l[1:NB])
    D[1:NB]=qq[1:NB]+2*dt*(SP[1:NB]+BP[1:NB]-DS[1:NB])

    # Top boundary conditions
    tau_ts = np.sqrt(tau_x**2 + tau_y**2) / rho
    A[0]=0
    B[0]=1.0
    C[0] = 0.0
    D[0] = (B_1**(2.0 / 3.0)) * tau_ts
    # Bottom boundary condition
    A[NB]=0.0;B[NB]=1.0;C[NB]=0;D[NB]=0
    tmp=solve_tri_implicit(A[0:NB], B[0:NB], C[0:NB], D[0:NB])
    qq_next[0:NB]=np.real(tmp)
    qq_next = np.maximum(qq_next, qq_min)
    if (NB < N):
        qq_next[NB+1:]=np.nan
#    a=qq_next-qq
#    print(a[::-1])
    return(qq_next)
