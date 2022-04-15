#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extracting time-differences and displaying results
"""
#%% Imports and functions

import time
t0 = time.time()
import numpy as np
from numpy import matlib


import supython as sup

import matplotlib as mpl
mpl.use('TkAgg')

from matplotlib import cm
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
#import matplotlib as mpl
import sys

from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

plt.rcParams['savefig.format'] = 'pdf'

#%% ------------------------- PLOTTING  FUNCTIONS ------------------------- %##

def show_subs(A,i=-1,dt=0.004,acaus=0,titles=None,suptitle=None,figsize=[6.4,4.8],length=5000,perc=100, intmethod = 'hanning',cmap=cm.seismic,name=None,subs=None,scaled=False,ylabel='Time [s]',xlabel='Lateral distance [m]',colorbar=None,eqclim=True,cmap2=cm.seismic,mask=None,shareax=1,clf=False):
    if isinstance(i, str) or (i > 0):
        fig=plt.figure(i,figsize=figsize)
        if clf:
            plt.clf()
    else:
        fig=plt.figure(figsize=figsize)
        
    if (name != None):
        fig.canvas.set_window_title(name)
    
    mats = len(A)
    if (subs==None):
        subs = (mats//10+1)*10 + (mats % 10)
    
    if (titles==None):
        titles=('',)*mats
    elif (len(titles) < mats):
        titles=titles + ('',)*(mats-len(titles))       
    
    if (mats > (subs//10)*(subs%10)):
        sys.exit('Amount of subplots smaller than amount of matrices')
    if (subs // 100 > 0):
        sys.exit('Maximum subplots size is 9 x 9')
        
    axes=np.array([None]*((subs//10)*(subs%10)))
    
    for k in np.arange(mats):        
        if np.all(np.isnan(A[k])):
            continue
        if (k == 0 or shareax == 0):
            axes[k] = fig.add_subplot(subs//10,subs%10,k+1)
        else:
            axes[k] = fig.add_subplot(subs//10,subs%10,k+1,sharex=axes[0],sharey=axes[0])
        if scaled:
            scale=np.nanmax(np.abs(A))
        else:
            scale=np.max(np.abs(A[k]))
        if (mask==None):
            im = show_image(A[k],dt=dt,acaus=acaus,i=None, intmethod =intmethod,title=titles[k],length=length,perc=perc,cmap=cmap,scale=scale,ylabel=ylabel,xlabel=xlabel,eqclim=eqclim) # minus to match gray scale of seismic unix
        else:
            im = show_image(A[k],mask=mask[k],cmap2=cmap2,dt=dt,acaus=acaus,i=None, intmethod =intmethod,title=titles[k],length=length,perc=perc,cmap=cmap,scale=scale,ylabel=ylabel,xlabel=xlabel,eqclim=eqclim)
    if (colorbar != None):
        fig.colorbar(im, ax=(axes[axes!=None]).tolist(), format=colorbar)
    if (suptitle != None): plt.suptitle(suptitle,fontsize=16,fontweight='bold')
        
def show_image(A,acaus=1,dt=0.004,i=-1,title='',ylabel='Time [s]',xlabel='Lateral distance [m]',intmethod = 'hanning',length=5000,perc=100,cmap=cm.gray,scale=None,eqclim=True,mask=None,cmap2=cm.gray):
    if (i == None):
        pass
    elif (i > 0):
        plt.figure(i)
        plt.clf()
    else:
        plt.figure()
    
    if not(scale==None) and not(eqclim):
        print('Only one of scale/eqlim can be used at a time, using matrix clim')
    
    if (scale == None or not(eqclim)):
        vmax=(perc/100)*abs(A).max()
        vmin=(perc/100)*abs(A).min()
    else:
        vmax=(perc/100)*scale
        vmin=-(perc/100)*scale
    
    if eqclim:
        if np.all(mask != None):
            plt.imshow(mask,cmap=cmap2,zorder=5)
        im = plt.imshow(A,cmap=cmap,vmax=vmax, vmin=-vmax,interpolation=intmethod)
        
    else:
        if np.all(mask != None):
            plt.imshow(mask,cmap=cmap2,zorder=5)
        im = plt.imshow(A,cmap=cmap,vmax=vmax, vmin=vmin,interpolation=intmethod)

    nt = A.shape[0]
    nx = A.shape[1]

    if acaus:
        tistart = np.floor(nt-np.floor(nt/100)*100)/2
        tiend = tistart + np.floor(nt/100)*100
        tdiff = np.floor((np.floor(nt/100)*100)/2)
        plt.yticks(np.linspace(tistart,tiend,11),np.round(np.linspace(-tdiff*dt,tdiff*dt,11),decimals=1))
        plt.plot([0,nx],[nt//2,nt//2],'k',linewidth=0.5)
    else:
        plt.yticks(np.linspace(0,np.floor(nt/100)*100,21),np.round(np.linspace(0,np.floor(nt/100)*100*dt,21),decimals=1))
    
    plt.xticks(np.linspace(0,nx-1,7),(np.round(np.linspace(-length/2,length/2,7),int(-np.floor(np.log(length))+2)),np.linspace((-length/2),(length/2),7,dtype=int))[length.is_integer() if isinstance(length,float) else True])
        

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
        
    plt.title(title)
    plt.axis('auto')
    plt.xlim([0,nx])

    return im
    
#%% Load Model

cp_b,hdr = sup.readsu("models/baseline_cp.su")
cp_m,hdr = sup.readsu("models/monitor_cp.su")
ro_b,hdr = sup.readsu("models/baseline_ro.su")
ro_m,hdr = sup.readsu("models/monitor_ro.su")

cp_sm,hdr = sup.readsu("models/smooth_cp.su")

nt = int(hdr.ns[0])
dz = np.round(hdr.d1[0],4)
ncell = hdr.ns.shape[0]
dxc = hdr.d2[0]

f2 = hdr.f2

diff = cp_m-cp_b

diff_res = diff.copy()
diff_res[diff < 100] = 0

cp = cp_b

outline = np.zeros_like(diff)
for i in range(diff.shape[0]-1):
    for j in range(diff.shape[1]-1):
        if diff[i,j] != diff[i+1,j]:
            outline[i,j] = 1
        if diff[i,j] != diff[i,j+1]:
            outline[i,j] = 1

color1 = colorConverter.to_rgba('white')

cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color1,color1],256)

cmap._init()
alphas = np.linspace(0,1,cmap.N+3)
cmap._lut[:,-1] = alphas

x_coords,y_coords = np.meshgrid(np.arange(diff.shape[1]),np.arange(diff.shape[0]))


#%% Get zero-offset times

p1ind = np.argmax((cp==3600)*np.tile(np.arange(641),[2401,1]).T,axis=0) + 10
mask1=np.zeros_like(cp)

for i in range(2401):
    mask1[:p1ind[i],i] = 1
  
p2ind = np.argmax((cp==2000)*np.tile(np.arange(641),[2401,1]).T,axis=0) + 20

mask2=np.zeros_like(cp)

for i in range(2401):
    mask2[:p2ind[i],i] = 1
    
p1time = np.sum((mask1)*2.5*2/cp,axis=0)[range(0,2401,4)]
p1time = gaussian_filter1d(p1time,sigma=25)
p2time = np.sum((mask2)*2.5*2/cp,axis=0)[range(0,2401,4)] 
p2time = gaussian_filter1d(p2time,sigma=25)
m1time = 2*p2time - 1*p1time 
m2time = 3*p2time - 2*p1time  
m3time = 4*p2time - 3*p1time  


#%%
def absmax(A):
    B = np.max(np.abs(A),axis=0)
    return B

def ricker(f0,nt,dt):
    t	=	np.arange(-nt/2*dt,(nt/2-1)*dt,dt)
    wavsym	= np.array(1-2*np.pi**2*f0**2*t**2)*np.exp(-np.pi**2*f0**2*t**2)
    wav	=	np.roll(wavsym,int(nt/2))
    return wav

def conv_wav(A,f0=25,dt=0.004):
    A_f = np.fft.rfft(A,axis=0)
    if f0 >0:
        wave = ricker(f0, A.shape[0], dt)
        wave_f = np.repeat(np.fft.rfft(wave,n=A.shape[0])[:,None],A.shape[1],axis=1)
    else:
        wave_f=1
    return np.fft.irfft(A_f*wave_f,axis=0)

def CrossCorr(A,B,subsamp=1,intp='spline'):
    if subsamp == 1:
        Af = np.fft.fft(A,axis=0)
        Bf = np.fft.fft(B,axis=0)
    elif intp=='spline':
        x = np.arange(A.shape[0])
       
        fa = interp1d(x, A, kind='cubic',axis=0)
        fb = interp1d(x, B, kind='cubic',axis=0)
        
        xnew = np.linspace(0,A.shape[0]-1,A.shape[0]*subsamp)
        
        Af = np.fft.fft(fa(xnew),axis=0)
        Bf = np.fft.fft(fb(xnew),axis=0) 
    elif intp=='freq':
        Af = np.zeros(( (A.shape[0]*subsamp,A.shape[1])),dtype=complex)
        Bf = np.zeros(( (B.shape[0]*subsamp,B.shape[1])),dtype=complex)
        Af[:int(A.shape[0]/2),:] = np.fft.fft(A,axis=0)[:int(A.shape[0]/2),:]
        Af[-int(A.shape[0]/2):,:] = np.fft.fft(A,axis=0)[-int(A.shape[0]/2):,:]
        
        Bf[:int(A.shape[0]/2),:] = np.fft.fft(B,axis=0)[:int(A.shape[0]/2),:]
        Bf[-int(A.shape[0]/2):,:] = np.fft.fft(B,axis=0)[-int(A.shape[0]/2):,:]
    else:
        raise('INVALID INTERPOLATION')
    
    CC = np.fft.fftshift(np.real(np.fft.ifft(Af*Bf.conj(),axis=0)),axes=(0))
    return CC/absmax(CC)

def R2I(A):
    return int(round(A))

#%% Load su files

b1,hdr = sup.readsu('Rb0z0_vmar_baseline.su')#
b2 = sup.readsuamp('Rbc0z0_vmar_baseline.su')#
b3 = sup.readsuamp('Rabc0z0_vmar_baseline.su')#

m1 = sup.readsuamp('Rb0z0_vmar_monitor.su')#
m2 = sup.readsuamp('Rbc0z0_vmar_monitor.su')#
m3 = sup.readsuamp('Rabc0z0_vmar_monitor.su')#

dtx=1
lnw=int(25)
m_alpha=2
samp=4
gaus=5
intp='freq'

f0 = 40
nt = int(hdr.ns[0])
dt = np.round(hdr.d1[0],4)
nrec = hdr.ns.shape[0]
dx = hdr.d2[0]

# add wavelet
b1 = conv_wav(b1[::,:],f0=f0,dt=dt)
m1 = conv_wav(m1[::,:],f0=f0,dt=dt)

b2 = conv_wav(b2[::,:],f0=f0,dt=dt)
m2 = conv_wav(m2[::,:],f0=f0,dt=dt)

b3 = conv_wav(b3[::,:],f0=f0,dt=dt)
m3 = conv_wav(m3[::,:],f0=f0,dt=dt)

tvec=np.arange(0,b1.shape[0])*dt       

if m_alpha==0:
    pmw=p1time
elif m_alpha==1:
    pmw=p2time
elif m_alpha==2:
    pmw=m1time
elif m_alpha==3:
    pmw=m2time
elif m_alpha==4:
    pmw=m3time

stw1=np.rint(p1time/dt).astype(int) - R2I(lnw/2) #1.8 #lef nums: 0.95, 1.25
edw1=stw1 + lnw  #1.95 #lef nums: 1.15, 1.45

stw2=np.rint(p2time/dt).astype(int) - R2I(lnw/2) #1.8 #lef nums: 0.95, 1.25
edw2=stw2 + lnw  #1.95 #lef nums: 1.15, 1.45

stw3=np.rint(m1time/dt).astype(int) - R2I(lnw/2) #1.8 #lef nums: 0.95, 1.25
edw3=stw3 + lnw 

stw4=np.rint(m2time/dt).astype(int) - R2I(lnw/2) #1.8 #lef nums: 0.95, 1.25
edw4=stw4 + lnw 

transparency=0.1

color4 = colorConverter.to_rgba('blue')

cmap4 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color4,color4],256)
cmap4._init()
alphas = np.linspace(0,transparency,cmap.N+3)
cmap4._lut[:,-1] = alphas

color2 = colorConverter.to_rgba('red')

cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color2,color2],256)

cmap2._init()
alphas = np.linspace(0,transparency,cmap.N+3)
cmap2._lut[:,-1] = alphas

color3 = colorConverter.to_rgba('white')

cmap3 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color3,color3],256)

cmap3._init()
alphas = np.linspace(0,transparency,cmap.N+3)
cmap3._lut[:,-1] = alphas

color5 = colorConverter.to_rgba('orange')

cmap5 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color5,color5],256)

cmap5._init()
alphas = np.linspace(0,transparency,cmap.N+3)
cmap5._lut[:,-1] = alphas

# fig = plt.figure(3)
ind=0
for arr in np.load('pickedwindows.npz').files:
    windows = ['edw1', 'edw2', 'edw3', 'edw4', 'stw1', 'stw2', 'stw3', 'stw4']
    globals()[windows[ind]] = np.load('pickedwindows.npz')[arr]
    ind+=1

ITER=0
def Plot_Windows(b1,num=10,interactive=True,lim1=1.6,lim2=0.):
    global edw1, edw2, edw3, edw4, stw1, stw2, stw3, stw4,dt
    x_coords,y_coords = np.meshgrid(np.arange(b1.shape[1]),np.arange(b1.shape[0]))
    if interactive:
        plt.figure("close figure after picking to continue")
    else:
        plt.figure(num)
    plt.clf()
    baselines=(b1,)

    i = 0
    if interactive:
        show_subs((baselines[i],),i="close figure after picking to continue",length=6000,cmap=cm.Greys,intmethod=None,eqclim=True,titles=('Baseline $R_{abc}^0$','Baseline $R^0_{bc}$','Baseline $R^0_b$'), perc=40,dt=dt,shareax=0)
    else:    
        show_subs((baselines[i],),i=num,length=6000,cmap=cm.Greys,intmethod=None,eqclim=True,titles=('Baseline $R_{abc}^0$','Baseline $R^0_{bc}$','Baseline $R^0_b$'), perc=40,dt=dt,shareax=0)
    
    ax = plt.gca()
   
    ax.xaxis.set_label_position('top'); ax.xaxis.tick_top();
    plt.xticks(np.linspace(0,nrec,7),np.arange(0,6001,1000),fontsize=14) 
    plt.xlabel('Lateral distance [$m$]', fontsize=18) 
    
    ax.set_title(('Baseline $R^0_{abc}$','Baseline $R^0_{bc}$','Baseline $R^0_b$')[ITER],fontsize=14)

    if interactive:
        plt.suptitle('Adjusting {:s} of {:s} {:d}'.format("top" if ind%2 == 0 else "bottom","primary" if ind//4 == 0 else "multiple",ind % 4 // 2 +1)) 
    ax.set_yticks(np.arange(0,1024,50))

    ax.set_yticklabels(np.round(np.arange(0,1001*dt,50*dt),decimals=1),fontsize=14)
    plt.ylabel('Time [$s$]', fontsize=18)
    
    plt.imshow(np.logical_and(y_coords>stw1,y_coords<edw1),cmap=cmap2,zorder=5)
    plt.imshow(np.logical_and(y_coords>stw2,y_coords<edw2),cmap=cmap3,zorder=5)
    plt.imshow(np.logical_and(y_coords>stw3,y_coords<edw3),cmap=cmap4,zorder=5)
    plt.imshow(np.logical_and(y_coords>stw4,y_coords<edw4),cmap=cmap5,zorder=5)
    plt.axis('auto')
    plt.ylim([lim1/dt,lim2/dt])
    plt.xlim([0,nrec])
    
    plt.contour(x_coords,y_coords,np.logical_and(y_coords>stw1,y_coords<edw1),colors='r',linewidths=.1)
    plt.contour(x_coords,y_coords,np.logical_and(y_coords>stw2,y_coords<edw2),colors='w',linewidths=.1)
    plt.contour(x_coords,y_coords,np.logical_and(y_coords>stw3,y_coords<edw3),colors='b',linewidths=.1)
    plt.contour(x_coords,y_coords,np.logical_and(y_coords>stw4,y_coords<edw4),colors='orange',linewidths=.1)
    plt.tight_layout()
    plt.draw()
    
plt.close('all')
fig = plt.figure("close figure after picking to continue")
ind=0

def on_close(event):
    global finishedpicking
    finishedpicking = True
    print('Windows are picked!')

def onclick(event):
    global ind, edw1, edw2, edw3, edw4, stw1, stw2, stw3, stw4
    
    windows = ['stw1', 'edw1', 'stw2', 'edw2', 'stw3', 'edw3', 'stw4', 'edw4']
    if event.button == 3:
        ind += 1
        if ind > 7:
            ind=0
    elif event.button == 1:
        newwindow = globals()[windows[ind]].copy()
        xm = int(event.xdata)
        xl = xm-25 if (xm-25)>0 else 0
        xr = xm+25 if (xm+25)<len(newwindow)-1 else len(newwindow)-1
        
        x = np.array([xl,xl+1, xm, xr-1,xr])
        y = np.array([newwindow[xl],newwindow[xl+1],int(event.ydata),newwindow[xr-1],newwindow[xr]])
        
        f = interp1d(x, y,kind='quadratic',axis=0)
        
        xnew = np.arange(xl,xr)
        ynew = f(xnew)
        
        newwindow[xl:xr] = ynew
        
        globals()[windows[ind]] = newwindow
        
    elif event.button ==2:
        try:
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                    ('double' if event.dblclick else 'single', event.button,
                    event.x, event.y, event.xdata, event.ydata))
        except:
            pass
        
    Plot_Windows(b1,lim1=1.6,lim2=0.4)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
fig.canvas.mpl_connect('close_event', on_close)
Plot_Windows(b1,lim1=1.6,lim2=0.4)
while not 'finishedpicking' in locals():
    plt.pause(5)

# np.savez('pickedwindows.npz',edw1, edw2, edw3, edw4, stw1, stw2, stw3, stw4)
del finishedpicking

Plot_Windows(b1,num=1,interactive=False)
plt.savefig('Figure5a') 
Plot_Windows(b2,num=2,interactive=False)
plt.savefig('Figure5b') 
Plot_Windows(b3,num=3,interactive=False)
plt.savefig('Figure5c') 

#%%

fig = plt.figure(4,figsize=[10.5,4.8])
#plt.savefig('Figure4.eps')


show_subs((cp_b,ro_b,),subs=13,i=4,dt=dz,length=12000,cmap=cm.viridis,intmethod=None,ylabel='Depth [m]',eqclim=False,scaled=True,shareax=0)
show_subs((np.nan,np.nan,diff,),subs=13,i=4,dt=dz,length=12000,intmethod=None,ylabel='Depth [m]',eqclim=True,perc=25)

for ax in fig.axes:
    ax.xaxis.set_label_position('top'); ax.xaxis.tick_top();
    ax.plot((0, 2401),(660/2.5, 660/2.5),'k--',lw=1)
    ax.plot((0, 2401),(1100/2.5, 1100/2.5),'k--',lw=1)
    
    
    
    plt.sca(ax)    
    ax.set_xticks(np.linspace(0,2401,4))
    ax.set_xticklabels(np.arange(0,6001,2000),fontsize=14)    
    plt.xlabel('Lateral distance [$m$]', fontsize=18) 
    
    cax = plt.colorbar(format='%d',orientation="horizontal", pad=0.01)
    if ax == fig.axes[1]:
        cax.set_label('Density $[kg/m^3]$   ',fontsize=18)
    else:
        cax.set_label('Velocity $[m/s]$  ',fontsize=18)
    
    if ax != fig.axes[2]:
        plt.contour(x_coords,y_coords,diff_res,colors='w',linewidths=.5) 
    else:
        cax.set_ticks([-25,0,25])
        cax.set_ticklabels([-25,0,100])
    cax.ax.tick_params(labelsize=14) 
    
    if ax == fig.axes[0]:
        ax.set_yticks(np.linspace(0,640,9))
        ax.set_yticklabels(np.arange(0,1601,200),fontsize=14)
        plt.ylabel('Depth [$m$]', fontsize=18)
    else:
        ax.set_yticks(np.linspace(0,640,9))
        ax.set_yticklabels(['']*9,fontsize=14)   
        plt.ylabel('')
plt.tight_layout()

plt.savefig('Figure4abc')         

thickness = np.sum(diff>0,axis=0)*dz

vel_b=cp_b[diff>0][0]
vel_m=cp_m[diff>0][0]

delt=(1/vel_b - 1/vel_m)*2*thickness*1000


#%%

ITER=1
for m_alpha in (1,2,3):
    if m_alpha==0:
        stw = stw1
        edw = edw1
    elif m_alpha==1:
        stw = stw2
        edw = edw2
    elif m_alpha==2:
        stw = stw3
        edw = edw3
    elif m_alpha==3:
        stw = stw4
        edw = edw4

    plt.figure(4+ITER,figsize=[9.6,3.6])#
    plt.clf()
    baselineresponses = (b3,b2,b1)
    monitorresponses = (m3,m2,m1)
    for J in range(3):
        blr = baselineresponses[J]
        mtr = monitorresponses[J]
        P1maskb=np.zeros_like(blr)
        Mmaskb=np.zeros_like(blr)
        
        P1maskm=np.zeros_like(blr)
        Mmaskm=np.zeros_like(blr)
        
        for i in range(nrec):
            P1maskb[stw1[i]:edw1[i],i] = blr[stw1[i]:edw1[i],i]
            Mmaskb[stw[i]:edw[i],i] = blr[stw[i]:edw[i],i]
            P1maskm[stw1[i]:edw1[i],i] = mtr[stw1[i]:edw1[i],i]
            Mmaskm[stw[i]:edw[i],i] = mtr[stw[i]:edw[i],i]   

        P1cut = P1maskb
        Mcut = Mmaskb
        
        CCb = CrossCorr(Mcut,P1cut,subsamp=samp,intp=intp)

        P1cut = P1maskm
        Mcut = Mmaskm
        
        CCm = CrossCorr(Mcut,P1cut,subsamp=samp,intp=intp)

        CC = CrossCorr(CCm,CCb,subsamp=1,intp=intp)
        
        if samp==1:
            tcaus = np.round(-((P1cut.shape[0]/2 if (P1cut.shape[0] % 2 == 0) else (P1cut.shape[0]-1)/2))*dt+np.arange(0,P1cut.shape[0]*samp)*dt/samp,3)
        else: 
            tcaus = np.round(-((P1cut.shape[0]/2 if (P1cut.shape[0] % 2 == 0) else (P1cut.shape[0])/2))*dt+np.arange(0,P1cut.shape[0]*samp)*dt/samp,3)
        
        Delta_twt=np.zeros(nrec)
        for i in range(nrec):
            Delta_twt[i] = tcaus[np.argmax(CC,axis=0)[i]]
            
        
        Delta_twt[np.abs(Delta_twt)*1000>15] = np.NaN
        
        if gaus:
            V=Delta_twt.copy()
            V[np.isnan(Delta_twt)]=0
            VV=gaussian_filter1d(V,sigma=gaus)
            
            W=np.ones_like(Delta_twt)
            W[np.isnan(Delta_twt)]=0
            WW=gaussian_filter1d(W,sigma=gaus)
            
            Delta_twt=VV/WW
            
        twt=np.zeros(nrec)
        for i in range(nrec):
            twt[i] = tcaus[np.argmax(CCb,axis=0)[i]]
            
        if gaus:
            V=twt.copy()
            V[np.isnan(twt)]=0
            VV=gaussian_filter1d(V,sigma=gaus)
            
            W=np.ones_like(twt)
            W[np.isnan(twt)]=0
            WW=gaussian_filter1d(W,sigma=gaus)
            
            twt=VV/WW
                
        refls = (r'$R_{abc}$',r'$R_{bc}$',r'$R_{b}$')
        subs = [np.NaN,]*9
        
        subs[J*3:J*3+2] = [CCb, CCm]
        show_subs(subs,subs=33,cmap=cm.gray_r,dt=0.001,acaus=1,i=7+ITER,length=6000,shareax=1,titles=('',refls[0],'','',refls[1],'','',refls[2],''))
        plt.ylim([3049,2048])
                
        subs = [np.NaN,]*9
        subs[J*3+2] = CC
        show_subs(subs,subs=33,cmap=cm.gray_r,dt=0.001,acaus=1,i=7+ITER,length=6000)
        plt.ylim([2024,2062])
        
        plt.suptitle('Correlations for {:s} {:d}'.format('multiple' if (m_alpha>1) else 'primary',m_alpha-1 if (m_alpha>1) else m_alpha+1),fontsize=13)
        plt.plot(range(nrec),2048+(Delta_twt)*1000,lw=2)
        
        
        mng = plt.get_current_fig_manager()
        try:
            mng.window.showMaximized()
        except AttributeError:
            #mng.frame.Maximize(True)
			# Some backends have no window (e.g., Agg)
            pass

        plt.tight_layout()
        
        plt.figure(4+ITER)
        plt.plot(range(nrec),(Delta_twt)*1000,lw=2)
         
    plt.bar(np.arange(0,nrec-1+.25,0.25),-m_alpha*delt,width=1,fill=True,color='cyan')
    plt.xticks(ticks=np.linspace(0,600,7),labels=np.arange(0,6001,1000),fontsize=14)
    plt.yticks(fontsize=14)    
    plt.xlim([0,600])
    plt.ylim([-15,15])
    plt.xlabel('Lateral distance [$m$]', fontsize=18)
    plt.ylabel(r'$\Delta$t [$ms$]', fontsize=18)
    plt.title('Change in time for {:s} {:d}'.format('multiple' if (m_alpha>1) else 'primary',m_alpha-1 if (m_alpha>1) else m_alpha+1),fontsize=18)
    plt.legend(('CC result: $R^0_{abc}$','CC result: $R^0_{bc}$','CC result: $R^0_{b}$','Actual'),fontsize=14,ncol=2,loc=1)
    plt.axhline(0,color='gray',lw=.5)
    
    ITER += 1
    
    plt.tight_layout()
	
    plt.savefig('Figure6{}'.format('abc'[m_alpha-1])) 
    

#%%
for ITER in range(1,4):
    fig = plt.figure(7+ITER)
    axes=fig.axes
    ax = axes[3]
    ax.set_ylim([3049,2048])
    ax = axes[6]
    ax.set_ylim([3049,2048])
    fig = plt.figure(7+ITER)
    for ax in fig.axes:
        plt.sca(ax)
        plt.xticks(ticks=np.linspace(0,600,7),labels=np.arange(0,6001,1000))
        
plt.pause(1)
for ITER in range(1,4):   
    plt.figure(7+ITER)
    plt.tight_layout()

plt.show()
