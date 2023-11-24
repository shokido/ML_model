import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
#fname_in="../Fortran/out_case_nn.txt"
out_form="X11"
#out_form="png"
fname_in="out_case_nn.txt"
f=open(fname_in,"r")
lines=f.readlines()
f.close()

iline=0
line=lines[iline].split()
ntime=int(line[0])
iline+=1
line=lines[iline].split()
z=[float(i) for i in line]
iline+=1
times=[]
temp=np.zeros((ntime,len(z)))
salt=np.zeros((ntime,len(z)))
u=np.zeros((ntime,len(z)))
v=np.zeros((ntime,len(z)))
for itime in range(0,ntime):
    line=lines[iline].split()
    iline+=1
    times.append(dt.datetime(int(line[0]),int(line[1]),int(line[2]),int(line[3]),int(line[4]),int(line[5])))

    line=lines[iline].split()
    iline+=1
    temp_tmp=[float(i) for i in line]
    line=lines[iline].split()
    iline+=1
    salt_tmp=[float(i) for i in line]
    line=lines[iline].split()
    iline+=1
    u_tmp=[float(i) for i in line]
    line=lines[iline].split()
    iline+=1
    v_tmp=[float(i) for i in line]

    temp[itime,:]=np.copy(np.asarray((temp_tmp)))
    salt[itime,:]=np.copy(np.asarray((salt_tmp)))
    u[itime,:]=np.copy(np.asarray((u_tmp)))
    v[itime,:]=np.copy(np.asarray((v_tmp)))

times=np.asarray(times)
z=np.asarray(z)
x,y=np.meshgrid(times,z)
clevs_t=np.arange(10,16,0.2)
clevs_s=np.arange(32,36.2,0.2)
clevs_uv=np.arange(-0.3,0.4,0.05)
xticks=times[0::5]
xticks=times[0::30]
xticks=times[0::60]
print(times)
ylim=[70,0]
fig = plt.figure(figsize=(12, 8))
plt.subplot(2,2,1)
plt.contourf(x,y,np.transpose(temp),cmap="jet",extend="both",levels=clevs_t)
plt.xticks(xticks)
plt.ylim(ylim)
plt.title("Potential temperature")
plt.xlabel("Time")
plt.ylabel("Depth")
plt.colorbar()

plt.subplot(2,2,2)
plt.contourf(x,y,np.transpose(salt),cmap="jet",extend="both",levels=clevs_s)
plt.xticks(xticks)
plt.ylim(ylim)
plt.title("Salinity")
plt.xlabel("Time")
plt.ylabel("Depth")
plt.colorbar()
plt.subplot(2,2,3)
plt.contourf(x,y,np.transpose(u),cmap="jet",extend="both",levels=clevs_uv)
plt.xticks(xticks)
plt.ylim(ylim)
plt.title("Zonal velocity")
plt.xlabel("Time")
plt.ylabel("Depth")
plt.colorbar()
plt.subplot(2,2,4)
plt.contourf(x,y,np.transpose(v),cmap="jet",extend="both",levels=clevs_uv)
plt.xticks(xticks)
plt.ylim(ylim)
plt.title("Meridional velocity")
plt.xlabel("Time")
plt.ylabel("Depth")
plt.colorbar()

# plt.subplot(2,1,2)
# plt.plot(temp[:,-1])
plt.tight_layout(rect=[0,0,1,0.96])
if (out_form=="png"):
    plt.savefig("sample_MLmodel.png")
else:
    plt.show()
