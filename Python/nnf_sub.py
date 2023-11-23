import numpy as np
from ml_util import *
# Closure constants
A_1 = 1.18;A_2 = 0.665;B_1 = 24.0;B_2 = 15.0
C_1 = 0.137;C_2 = 0.75;C_3 = 0.352;C_5 = 0.2
SQ = 3.0
q2_min = 1.0e-8
q2l_min = 1.0e-8
gamma1 = 1.0 / 3.0 - 2.0 * A_1 / B_1
gamma2 = (2.0 * A_1 * (3.0 - 2.0 * C_2) + B_2 * (1.0 - C_3)) / B_1
F_1 = B_1 * (gamma1 - C_1) + 2.0 * A_1 * (3.0 - 2.0 * C_2) + 3.0 * A_2 * (1.0 - C_2) * (1.0 - C_5)
F_2 = B_1 * (gamma1 + gamma2) - 3.0 * A_1 * (1.0 - C_2)
R_f1 = B_1 * (gamma1 - C_1) / F_1
R_f2 = B_1 * gamma1 / F_2
R_fc = gamma1 / (gamma1 + gamma2)
tiny = 1.0e-12

def cal_Rf(Ri):
    R_i1 = 0.5 * A_2 * F_2 / (A_1 * F_1)
    R_i2 = 0.5 * R_f1 / R_i1
    R_i3 = (2.0 * R_f2 - R_f1) / R_i1
    Rf = R_i1 * (Ri + R_i2 - np.sqrt(Ri**2 - R_i3 * Ri + R_i2**2))
    Rf[Rf >= R_fc] = R_fc - tiny
    return Rf

def cal_SH_my2(Rf):
    return 3.0 * A_2 * (gamma1 + gamma2) * (R_fc - Rf) / (1.0 - Rf)

def cal_SM_my2(Rf):
    SH2 = cal_SH_my2(Rf)
    return (A_1 * F_1) / (A_2 * F_2) * ((R_f1 - Rf) / (R_f2 - Rf)) * SH2

def cal_SH_mynnf25(GH, GM, alpha_C):
    phi_5 = 1.0 - 3.0 * (alpha_C**2) * A_2 * B_2 * (1.0 - C_3) * GH
    phi_1 = 1.0 - 9.0 * (alpha_C**2) * A_1 * A_2 * (1.0 - C_2) * GH
    phi_2 = phi_5 + 9.0 * (alpha_C**2) * (A_2**2) * (1.0 - C_2) * (1.0 - C_5) * GH
    phi_3 = phi_5 - 12.0 * (alpha_C**2) * (A_1 * A_2) * (1.0 - C_2) * GH
    phi_4 = 6.0 * (alpha_C**2) * (A_1**2) * GM
    numer = alpha_C * A_2 * (phi_1 + 3.0 * C_1 * phi_4)
    denom = phi_1 * phi_3 + phi_2 * phi_4
    return numer / denom

def cal_SM_mynnf25(GH, GM, alpha_C):
    phi_5 = 1.0 - 3.0 * (alpha_C**2) * A_2 * B_2 * (1.0 - C_3) * GH
    phi_1 = 1.0 - 9.0 * (alpha_C**2) * A_1 * A_2 * (1.0 - C_2) * GH
    phi_2 = phi_5 + 9.0 * (alpha_C**2) * (A_2**2) * (1.0 - C_2) * (1.0 - C_5) * GH
    phi_3 = phi_5 - 12.0 * (alpha_C**2) * (A_1 * A_2) * (1.0 - C_2) * GH
    phi_4 = 6.0 * (alpha_C**2) * (A_1**2) * GM
    numer = alpha_C * A_1 * (phi_2 - 3.0 * C_1 * phi_3)
    denom = phi_1 * phi_3 + phi_2 * phi_4
    return numer / denom

def mynnf25_vdiff(N, bvf, shear2, Q2, L):
    GH = -1*bvf * L[1:N]**2 / Q2[1:N]
    GM = shear2 * L[1:N]**2 / Q2[1:N]
    Ri = bvf / shear2
    Rf = cal_Rf(Ri)
    SM2 = cal_SM_my2(Rf)
    alpha_C = np.sqrt(1.0 / np.maximum(tiny, (B_1 * cal_SM_my2(Rf) * (1.0 - Rf) * GM)))
    alpha_C = np.minimum(alpha_C, 1.0)
    SM = cal_SM_mynnf25(GH, GM, alpha_C)
    SH = cal_SH_mynnf25(GH, GM, alpha_C)
    KM = L[1:N] * np.sqrt(Q2[1:N]) * SM
    KT = L[1:N] * np.sqrt(Q2[1:N]) * SH
    KS = L[1:N] * np.sqrt(Q2[1:N]) * SH
    return KM, KT, KS


# Subroutine diag_l_mynnf
def diag_l_mynnf(N, Z_Q, bvf, Q2, tau_x, tau_y, sw, nsw, salflux):
    L = np.zeros(N+1)
    bvf_targ = np.zeros(N)
    mo_inv = cal_mo_inv(T0, S0, tau_x, tau_y, sw, nsw, salflux)
    L[N] = 0.0
    q_int = 0.0
    qz_int = 0.0
    dz = np.diff(Z_Q)
    func1 = np.sqrt(Q2)
    q_int = np.sum(dz * (func1[1:] + func1[:-1]) * 0.5)
    func1 = np.sqrt(Q2) * np.abs(Z_Q)
    qz_int = np.sum(dz * (func1[1:] + func1[:-1]) * 0.5)
    lt_inv = (q_int / (0.23 * qz_int)) * np.ones(N)
    xi = np.abs(Z_Q[0:N]) * mo_inv
    ls_inv=np.copy(xi)
    ind_1=np.where(xi>=1)[0]
    if (len(ind_1)>0):
        ls_inv[ind_1]=3.7/(kappa * np.abs(Z_Q[ind_1]))
    ind_2=np.where((xi>=0)&(xi<1))[0]
    if (len(ind_2)>0):
        ls_inv[ind_2]=(1.0 + 2.7 * xi[ind_2])/(kappa * np.abs(Z_Q[ind_2]))
    ind_3=np.where(xi<0)[0]
    if (len(ind_3)>0):
        ls_inv[ind_3]=((1.0 - 2.7 * xi[ind_2])**(-0.2))/(kappa * np.abs(Z_Q[ind_3]))
    bvf_targ[1:N]=np.copy(bvf)
    bvf_targ[0] = bvf[0]
    mask = bvf_targ >= 0
    lb_inv = np.where(mask, np.sqrt(bvf_targ / (Q2[0:N]+tiny*np.ones(N))) / 0.53, 0.0 + tiny)

    l_inv = ls_inv + lt_inv + lb_inv
    #print(ls_inv)
    L[0:N] = 1.0 / l_inv
    return L
# Subroutine cal_kq
def cal_kq(n, km):
    kq=np.zeros(n)
    kq[0] = SQ * km[0]
    kq[1:(n-1)]=0.5*(km[0:(n-2)]+km[1:(n-1)])*SQ
    kq[n-1]=0.5*km[n-2]*SQ
    return(kq)
