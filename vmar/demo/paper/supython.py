"""
- SUpython -
Scripts to read in work on and write out su files in Python
Created in July 2019 at Delft University of Technology
@author: Joeri Brackenhoff (J.A.Brackenhoff@tudelft.nl)
Additional contributions by Johno van IJsseldijk (J.E.vanIJsseldijk@tudelft.nl)
"""

# Import required modules
import struct
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from time import gmtime

# Define required global variables for all attributes of SU header and their type
vari = ['tracl','tracr','fldr','tracf','ep','cdp','cdpt','trid','nvs','nhs','duse','offset','gelev','selev',
        'sdepth','gdel','sdel','swdep','gwdep','scalel','scalco','sx','sy','gx','gy','counit','wevel','swevel',
        'sut','gut','sstat','gstat','tstat','laga','lagb','delrt','muts','mute','ns','dt','gain','igc','igi',
        'corr','sfs','sfe','slen','styp','stas','stae','tatyp','afilf','afils','nofilf','nofils','lcf','hcf',
        'lcs','hcs','year','day','hour','minute','sec','timbas','trwf','grnors','grnofr','grnlof','gaps',
        'otrav','d1','f1','d2','f2','ungpow','unscale','ntr','mark','shortpad',];
vari_bytes = 'iiiiiiihhhhiiiiiiiihhiiiihhhhhhhhhhhhhHHhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhffffffihhhhhhhhhhhhhhhh';
  
# Define all of the attributes of a SU header  
class suhdr:
    def __init__(self,hdr):
        self.tracl=hdr[0,:];
        self.tracr=hdr[1,:];
        self.fldr=hdr[2,:];
        self.tracf=hdr[3,:];
        self.ep=hdr[4,:];
        self.cdp=hdr[5,:];
        self.cdpt=hdr[6,:];
        self.trid=hdr[7,:];
        self.nvs=hdr[8,:];
        self.nhs=hdr[9,:];
        self.duse=hdr[10,:];
        self.offset=hdr[11,:];
        self.gelev=hdr[12,:];
        self.selev=hdr[13,:];
        self.sdepth=hdr[14,:];
        self.gdel=hdr[15,:];
        self.sdel=hdr[16,:];
        self.swdep=hdr[17,:];
        self.gwdep=hdr[18,:];
        self.scalel=hdr[19,:];
        self.scalco=hdr[20,:];
        self.sx=hdr[21,:];
        self.sy=hdr[22,:];
        self.gx=hdr[23,:];
        self.gy=hdr[24,:];
        self.counit=hdr[25,:];
        self.wevel=hdr[26,:];
        self.swevel=hdr[27,:];
        self.sut=hdr[28,:];
        self.gut=hdr[29,:];
        self.sstat=hdr[30,:];
        self.gstat=hdr[31,:];
        self.tstat=hdr[32,:];
        self.laga=hdr[33,:];
        self.lagb=hdr[34,:];
        self.delrt=hdr[35,:];
        self.muts=hdr[36,:];
        self.mute=hdr[37,:];
        self.ns=hdr[38,:];
        self.dt=hdr[39,:];
        self.gain=hdr[40,:];
        self.igc=hdr[41,:];
        self.igi=hdr[42,:];
        self.corr=hdr[43,:];
        self.sfs=hdr[44,:];
        self.sfe=hdr[45,:];
        self.slen=hdr[46,:];
        self.styp=hdr[47,:];
        self.stas=hdr[48,:];
        self.stae=hdr[49,:];
        self.tatyp=hdr[50,:];
        self.afilf=hdr[51,:];
        self.afils=hdr[52,:];
        self.nofilf=hdr[53,:];
        self.nofils=hdr[54,:];
        self.lcf=hdr[55,:];
        self.hcf=hdr[56,:];
        self.lcs=hdr[57,:];
        self.hcs=hdr[58,:];
        self.year=hdr[59,:];
        self.day=hdr[60,:];
        self.hour=hdr[61,:];
        self.minute=hdr[62,:];
        self.sec=hdr[63,:];
        self.timbas=hdr[64,:];
        self.trwf=hdr[65,:];
        self.grnors=hdr[66,:];
        self.grnofr=hdr[67,:];
        self.grnlof=hdr[68,:];
        self.gaps=hdr[69,:];
        self.otrav=hdr[70,:];
        self.d1=hdr[71,:];
        self.f1=hdr[72,:];
        self.d2=hdr[73,:];
        self.f2=hdr[74,:];
        self.ungpow=hdr[75,:];
        self.unscale=hdr[76,:];
        self.ntr=hdr[77,:];
        self.mark=hdr[78,:];
        self.shortpad=hdr[79,:];
        
# Define the attributes that are scaled to the actual data
class suhdrscale:
    def __init__(self,hdr):
        scalel=hdr[19,0];
        scalco=hdr[20,0];
        
        if scalel<0:
            scaledep=float(-1.0/scalel)
        elif scalel==0:
            scaledep=1.0;
        else:
            scaledep=float(scalel)
        
        if scalco<0:
            scalepos=float(-1.0/scalco)
        elif scalco==0:
            scalepos=1.0;
        else:
            scalepos=float(scalco)
        
        self.tracl=hdr[0,:];
        self.tracr=hdr[1,:];
        self.fldr=hdr[2,:];
        self.tracf=hdr[3,:];
        self.ep=hdr[4,:];
        self.cdp=hdr[5,:];
        self.cdpt=hdr[6,:];
        self.trid=hdr[7,:];
        self.nvs=hdr[8,:];
        self.nhs=hdr[9,:];
        self.duse=hdr[10,:];
        self.offset=hdr[11,:];
        self.gelev=np.array(hdr[12,:],dtype=float)*scaledep;
        self.selev=np.array(hdr[13,:],dtype=float)*scaledep;
        self.sdepth=np.array(hdr[14,:],dtype=float)*scaledep;
        self.gdel=np.array(hdr[15,:],dtype=float)*scaledep;
        self.sdel=np.array(hdr[16,:],dtype=float)*scaledep;
        self.swdep=np.array(hdr[17,:],dtype=float)*scaledep;
        self.gwdep=np.array(hdr[18,:],dtype=float)*scaledep;
        self.scalel=hdr[19,:];
        self.scalco=hdr[20,:];
        self.sx=np.array(hdr[21,:],dtype=float)*scalepos;
        self.sy=np.array(hdr[22,:],dtype=float)*scalepos;
        self.gx=np.array(hdr[23,:],dtype=float)*scalepos;
        self.gy=np.array(hdr[24,:],dtype=float)*scalepos;
        self.counit=hdr[25,:];
        self.wevel=hdr[26,:];
        self.swevel=hdr[27,:];
        self.sut=hdr[28,:];
        self.gut=hdr[29,:];
        self.sstat=hdr[30,:];
        self.gstat=hdr[31,:];
        self.tstat=hdr[32,:];
        self.laga=hdr[33,:];
        self.lagb=hdr[34,:];
        self.delrt=hdr[35,:];
        self.muts=hdr[36,:];
        self.mute=hdr[37,:];
        self.ns=hdr[38,:];
        self.dt=hdr[39,:];
        self.gain=hdr[40,:];
        self.igc=hdr[41,:];
        self.igi=hdr[42,:];
        self.corr=hdr[43,:];
        self.sfs=hdr[44,:];
        self.sfe=hdr[45,:];
        self.slen=hdr[46,:];
        self.styp=hdr[47,:];
        self.stas=hdr[48,:];
        self.stae=hdr[49,:];
        self.tatyp=hdr[50,:];
        self.afilf=hdr[51,:];
        self.afils=hdr[52,:];
        self.nofilf=hdr[53,:];
        self.nofils=hdr[54,:];
        self.lcf=hdr[55,:];
        self.hcf=hdr[56,:];
        self.lcs=hdr[57,:];
        self.hcs=hdr[58,:];
        self.year=hdr[59,:];
        self.day=hdr[60,:];
        self.hour=hdr[61,:];
        self.minute=hdr[62,:];
        self.sec=hdr[63,:];
        self.timbas=hdr[64,:];
        self.trwf=hdr[65,:];
        self.grnors=hdr[66,:];
        self.grnofr=hdr[67,:];
        self.grnlof=hdr[68,:];
        self.gaps=hdr[69,:];
        self.otrav=hdr[70,:];
        self.d1=hdr[71,:];
        self.f1=hdr[72,:];
        self.d2=hdr[73,:];
        self.f2=hdr[74,:];
        self.ungpow=hdr[75,:];
        self.unscale=hdr[76,:];
        self.ntr=hdr[77,:];
        self.mark=hdr[78,:];
        self.shortpad=hdr[79,:];
        self.t=np.linspace(self.f1[0],self.f1[0]+(self.ns[0]-1)*self.d1[0],int(self.ns[0]));

# Read in the header data and the amplitudes of the traces
def readsu(filename, scale=1):
    
    # filename = path to the file to be read in
    # scale    = scale the header data from SU format to python format (=1) or not (=0)
    
    # Set global variables
    global vari
    global vari_bytes
    
    # Check whether the file path exists
    filetest = Path(filename)
    if filetest.is_file()==False:
        raise Exception('File could not be found')    
    
    # Open the data to determine the size
    fp = open(filename,'rb');
    fp.seek(0,2)
    size = fp.tell()
    
    # Read in the first header
    fp.seek(0,0)
    file = fp.read(240);
    inivalues = struct.unpack(vari_bytes,file)
    
    # From the header grab the trace samples and the amount of receivers and allocate the data
    ns = inivalues[38]
    nx = int(size/(240+4*ns))
    amp = np.zeros((ns,nx))
    
    # Allocate the header object
    nrval = len(vari)
    values = np.zeros((nrval,nx))
    values[:,0] = inivalues[0:nrval]
    
    # Set the format for reading in one trace
    fmt = "<%df" % ns;
    
    # Loop over the file to read in the traces and the headers
    for ix in range(nx):
        file = fp.read(4*ns)
        tmp = np.reshape(struct.unpack(fmt,file),(ns))
        amp[:,ix] = tmp;
        file = fp.read(240);
        if ix<(nx-1):
            inivalues = struct.unpack(vari_bytes,file)
            values[:,ix+1] = inivalues[0:nrval]
           
    # Close the file
    fp.close()
    
    # Convert the header data to the hdr class and scale the data if asked
    if scale==0:
        hdr = suhdr(values);
    elif scale==1:
        hdr = suhdrscale(values)
        
    return(amp,hdr)
    
# Read in the amplitudes of the traces
def readsuamp(filename):        
    
    # filename = path to the file to be read in
    
    # Set global variables
    global vari_bytes
    
    # Check whether the file path exists
    filetest = Path(filename)
    if filetest.is_file()==False:
        raise Exception('File could not be found')
    
    # Open the data to determine the size
    fp = open(filename,'rb');
    fp.seek(0,2)
    size = fp.tell()
    
    # Read in the first header
    fp.seek(0,0)
    file = fp.read(240);
    inivalues = struct.unpack(vari_bytes,file)
    
    # From the header grab the trace samples and the amount of receivers and allocate the data
    ns = inivalues[38]
    nx = int(size/(240+4*ns))
    amp = np.zeros((ns,nx))
    
    # Set the format for reading in one trace
    fmt = "<%df" % ns;
    
    # Loop over the file to read in the traces
    for ix in range(nx):
        file = fp.read(4*ns)
        tmp = np.reshape(struct.unpack(fmt,file),(ns))
        amp[:,ix] = tmp;
        fp.seek(240,1)
           
    # Close the file
    fp.close()
        
    return(amp)
    
# Read in the header data
def readsuhdr(filename, scale=1):
    
    # filename = path to the file to be read in
    # scale    = scale the header data from SU format to python format (=1) or not (=0)
    
    # Set global variables
    global vari
    global vari_bytes
    
    # Check whether the file path exists
    filetest = Path(filename)
    if filetest.is_file()==False:
        raise Exception('File could not be found')
    
    
    # Open the data to determine the size
    fp = open(filename,'rb');
    fp.seek(0,2)
    size = fp.tell()
    
    # Read in the first header
    fp.seek(0,0)
    file = fp.read(240);
    inivalues = struct.unpack(vari_bytes,file)
    
    # From the header grab the trace samples and the amount of receivers
    ns = inivalues[38]
    nx = int(size/(240+4*ns))
    
    # Allocate the header object
    nrval = len(vari)
    values = np.zeros((nrval,nx))
    values[:,0] = inivalues[0:nrval]
    
    # Loop over the file to read in the headers
    for ix in range(nx):
        fp.seek(ns*4,1)
        file = fp.read(240);
        if ix<(nx-1):
            inivalues = struct.unpack(vari_bytes,file)
            values[:,ix+1] = inivalues[0:nrval]
           
    # Close the file
    fp.close()
    
    # Convert the header data to the hdr class and scale the data if asked
    if scale==0:
        hdr = suhdr(values);
    elif scale==1:
        hdr = suhdrscale(values)
        
    return(hdr)
    
# Write out data from python format to SU format
def writesu(filename, amp, hdr, scale=1):
    
    # filename = path to the file to be written out
    # amp      = matrix containing the amplitudes of the data
    # hdr      = header object containing the headers for the data
    # scale    = scale the header data from python format to SU format (=1) or not (=0)
    
    # Set global variables
    global vari
    global vari_bytes
    
    # Scale the data back to the SU format
    if scale==1:
        scalel=hdr.scalel[0];
        scalco=hdr.scalco[0];
        
        if scalel<0:
            scaledep=float(-1.0/scalel)
        elif scalel==0:
            scaledep=1.0;
        else:
            scaledep=float(scalel)
        
        if scalco<0:
            scalepos=float(-1.0/scalco)
        elif scalco==0:
            scalepos=1.0;
        else:
            scalepos=float(scalco)
    else:
        scaledep=1.0;
        scalepos=1.0;
    
    # Open the file to write data and get the size of the data to be written
    fpout = open(filename,'wb');
    [nz,nx] = np.shape(amp)

    # Write out the data, check for the appropiate header format for every attribute
    for ix in range(nx):
        for ih in range(94):
            if ih < 80:
                if vari[ih] in vari[12:19]:
                    attrib = getattr(hdr,vari[ih])/scaledep;
                elif vari[ih] in vari[21:25]:
                    attrib = getattr(hdr,vari[ih])/scalepos;
                else:
                    attrib = getattr(hdr,vari[ih])
            if vari_bytes[ih] == 'i':
                fpout.write(struct.pack(vari_bytes[ih],int(attrib[ix])))
            elif vari_bytes[ih] == 'h' and ih < 80:
                fpout.write(struct.pack(vari_bytes[ih],np.short(attrib[ix])))
            elif vari_bytes[ih] == 'h' and ih > 79:
                fpout.write(struct.pack(vari_bytes[ih],np.short(0)))
            elif vari_bytes[ih] == 'H':
                fpout.write(struct.pack(vari_bytes[ih],np.ushort(attrib[ix])))
            elif vari_bytes[ih] == 'f':
                fpout.write(struct.pack(vari_bytes[ih],float(attrib[ix])))
        data = amp[:,ix]
        fpout.write(struct.pack('<%df' % len(data), *data))
        
    # Close the data
    fpout.close()
 
# Create a header object and set the most important values in the SU format
def makehdr(amp, dx=10, dt=0.004, t0=0, f2=-3000, scl=-1000, gelev=0, sdepth=0):
    
    # Determine the amount of samples, the amount of receivers and the amount of shots
    ns = amp.shape[0]
    nx = amp.shape[1]
    shots = 1 if len(amp.shape) == 2 else amp.shape[2]
    
    # Allocate the header size
    outsize = np.ones((nx*shots))
    
    # Create the header object
    hdr = suhdr
    
    # Set the header values
    hdr.tracl=np.tile(np.arange(1,1+nx),shots);
    hdr.tracr=outsize*0;
    hdr.fldr=outsize if shots==1 else np.repeat(np.arange(1,1+shots),nx);
    hdr.tracf=np.arange(1,1+nx*shots);
    hdr.ep=outsize*0;
    hdr.cdp=outsize*0;
    hdr.cdpt=outsize*0;
    hdr.trid=outsize;
    hdr.nvs=outsize*0;
    hdr.nhs=outsize*0;
    hdr.duse=outsize*0;
    hdr.d1=outsize*dt;
    hdr.f1=outsize*t0;
    hdr.d2=outsize*dx;
    hdr.f2=outsize*f2;
    hdr.gelev=outsize*gelev;
    hdr.sdepth=outsize*sdepth;
    hdr.selev=-hdr.sdepth;
    hdr.gdel=outsize*0;
    hdr.sdel=outsize*0;
    hdr.swdep=outsize*0;
    hdr.gwdep=outsize*0;
    hdr.scalel=outsize*scl;
    hdr.scalco=outsize*scl;
    hdr.sx=outsize*0 if shots==1 else np.repeat(np.arange(0,shots*dx,dx),nx)+f2;
    hdr.sy=outsize*0;
    hdr.gx=np.tile(np.arange(0,nx*dx,dx),shots)+f2;
    hdr.gy=outsize*0;
    hdr.offset=hdr.gx-hdr.sx;
    hdr.counit=outsize*0;
    hdr.wevel=outsize*0;
    hdr.swevel=outsize*0;
    hdr.sut=outsize*0;
    hdr.gut=outsize*0;
    hdr.sstat=outsize*0;
    hdr.gstat=outsize*0;
    hdr.tstat=outsize*0;
    hdr.laga=outsize*0;
    hdr.lagb=outsize*0;
    hdr.delrt=outsize*0;
    hdr.muts=outsize*0;
    hdr.mute=outsize*0;
    hdr.ns=outsize*ns;
    hdr.dt=outsize*dt*1e6;
    hdr.gain=outsize*0;
    hdr.igc=outsize*0;
    hdr.igi=outsize*0;
    hdr.corr=outsize*0;
    hdr.sfs=outsize*0;
    hdr.sfe=outsize*0;
    hdr.slen=outsize*0;
    hdr.styp=outsize*0;
    hdr.stas=outsize*0;
    hdr.stae=outsize*0;
    hdr.tatyp=outsize*0;
    hdr.afilf=outsize*0;
    hdr.afils=outsize*0;
    hdr.nofilf=outsize*0;
    hdr.nofils=outsize*0;
    hdr.lcf=outsize*0;
    hdr.hcf=outsize*0;
    hdr.lcs=outsize*0;
    hdr.hcs=outsize*0;
    hdr.year=outsize*gmtime().tm_year;
    hdr.day=outsize*gmtime().tm_yday;
    hdr.hour=outsize*gmtime().tm_hour;
    hdr.minute=outsize*gmtime().tm_min;
    hdr.sec=outsize*gmtime().tm_sec;
    hdr.timbas=outsize*4;
    hdr.trwf=outsize*nx;
    hdr.grnors=outsize*0;
    hdr.grnofr=outsize*0;
    hdr.grnlof=outsize*0;
    hdr.gaps=outsize*0;
    hdr.otrav=outsize*0;
    hdr.ungpow=outsize*0;
    hdr.unscale=outsize*0;
    hdr.ntr=outsize*shots*nx;
    hdr.mark=outsize*0;
    hdr.shortpad=outsize*0;
    hdr.t=np.linspace(hdr.f1[0],hdr.f1[0]+(hdr.ns[0]-1)*hdr.d1[0],int(hdr.ns[0]));
    
    return hdr
    
# Plot data using the SU hdr info
def plotsu(Z,X,Y,comap='gray',vmin=0.0,vmax=0.0,xlabel='',ylabel='',clabel='',asp=1,cm=1):
    
    # Z      = Amplitudes of the data
    # X      = vector containing the x-values of the extent of the data
    # Y      = vector containing the y-values of the extent of the data
    # comap  = Colormap for the plotting, standard is grayscale
    # vmin   = minimum for the colorbar
    # vmax   = maximum for the colorbar
    # xlabel = label for the x-axis
    # ylabel = label for the y-axis
    # clabel = label for the colorbar
    # asp    = aspect ratio
    # cm     = colorbar on (=1) or off (=0)
    
    # Set the extent for the image
    extent=[X[0],X[-1],Y[-1],Y[0]]
    
    # If the minimum and maximum values are not chosen, pick them from the data
    if (vmin==0.0):
        vmin1 = np.min(Z)
    else:
        vmin1 = vmin
    if (vmax==0.0):
        vmax1 = np.max(Z)
    else:
        vmax1 = vmax
    
    # Plot the data
    plt.imshow(Z,cmap=comap,extent=extent,aspect=asp,vmin=vmin1,vmax=vmax1)
    # Plot colorbar if set
    if cm==1:
        cbar=plt.colorbar()
        cbar.set_label(clabel, rotation=270,labelpad=20)
    # Set extent and labels
    plt.xlim([extent[0],extent[1]])
    plt.ylim([extent[3],extent[2]])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.gca().invert_yaxis()