# -*- coding: utf-8 -*-
"""
Created on Sat May 23 12:26:06 2015

@author: xhxiao
"""

#from tomopy.io.reader import *
#import numpy as np
#import os.path
#from numpy.testing import assert_allclose
import tomopy 
from tomopy.recon.rotation import write_center
from tomopy.recon.algorithm import recon
from tomopy import minus_log
import os, h5py, glob, fnmatch
import numpy as np
from scipy import misc
import time

####----------------- Section 1: Input file -- Staragt ----------------- 
print time.asctime() 
data_top_dir = '/media/2BM_Backup46/2016_06/Commissioning/'
dir_idx = 43

data_dir = glob.glob(os.path.join(data_top_dir,'Exp'+'{:03}'.format(dir_idx)+'*'))[0]
output_dir = data_dir

data_files = np.sort(fnmatch.filter(os.listdir(data_dir),'*.hdf'))

file_data = data_files[0]
file_flat = data_files[1]
file_dark = data_files[2]
print file_data
print file_flat
print file_dark

file_name = os.path.join(data_dir, file_data)
flat_name = os.path.join(data_dir, file_flat)
dark_name = os.path.join(data_dir, file_dark)
output_file = output_dir+'/recon_'+file_data.split(".")[-2]+'/recon_'+file_data.split(".")[-2] +'_'
####----------------- Section 1: Input file -- End -----------------


####----------------- Section 2: Parameter configuration -- Start ----------------- 
f = h5py.File(file_name,"r")
try:
    arr = f["exchange/data"]
except:
    try:
        arr = f["entry/data/data"]
    except:
        arr = f["MyDataDir/TomoShots"] 
numProj = (arr.shape)[0]
numSlices = (arr.shape)[1]
imgWidth = (arr.shape)[2]
f.close()

offset = 0
numRecSlices = numSlices
chunk_size = 300
if chunk_size > numRecSlices:
    chunk_size = numRecSlices
margin_slices = 30
num_chunk = np.int(numRecSlices/(chunk_size-margin_slices)) + 1
if numRecSlices == chunk_size:
    num_chunk = 1
  
z = 6
eng = 27.4
pxl = 0.65e-4
rat = 20e-04    

zinger_level = 500
####----------------- Section 2: Parameter configuration -- End -----------------   


########----------------- Section 3: Finding center -- Start -----------------
#sliceStart = 300
#sliceEnd = 320
#f = h5py.File(file_name,"r")
#try:
#    arr = f["exchange/data"]
#except:
#    try:
#        arr = f["entry/data/data"]
#    except:
#        arr = f["MyDataDir/TomoShots"]     
#data = arr[1:numProj,sliceStart:sliceEnd,:]
#f.close()
#
#f = h5py.File(flat_name,"r")
#try:
#    arr = f["exchange/data_white"]
#except:
#    try:
#        arr = f["entry/data/data"]
#    except:
#        arr = f["MyDataDir/TomoShots"] 
#white = arr[1:9,sliceStart:sliceEnd,:]
#f.close()
#
#f = h5py.File(dark_name,"r")
#try:
#    arr = f["exchange/data_dark"]
#except:
#    try:
#        arr = f["entry/data/data"]
#    except:
#        arr = f["MyDataDir/TomoShots"]     
#dark = arr[1:9,sliceStart:sliceEnd,:]
#f.close()
#
#data_size = data.shape
#theta = np.linspace(0,np.pi,num=data_size[0]+1)
#
## zinger_removal
#data = tomopy.misc.corr.remove_outlier(data,zinger_level,size=15,axis=0)
#white = tomopy.misc.corr.remove_outlier(white,zinger_level,size=15,axis=0)
#                           
## there is modification in below normalize routine
#data = tomopy.prep.normalize.normalize(data,white,dark)
#data = tomopy.prep.normalize.normalize_bg(data,air=10)
#
## strip removal
#tomopy.set_debug('True')
#data = tomopy.prep.stripe.remove_stripe_fw(data,level=10,wname='sym16',sigma=1,pad=True)
#data = tomopy.prep.phase.retrieve_phase(data,pixel_size=pxl,dist=z,energy=eng,alpha=rat,pad=True)
#tomopy.set_debug('False')
    
#data_size = data.shape
# 
#center_shift = 0
#center_shift_w = 10
## there is modification in the below routine
#data = tomopy.prep.normalize.minus_log(data)
#write_center(data[:,9:11,:], theta, dpath=output_dir+'/data_center/', 
#             cen_range=(data.shape[2]/2+center_shift,data.shape[2]/2+center_shift+center_shift_w,0.5),emission=False)
########----------------- Section 3: Finding center -- End -----------------

center = 1301.5		#Exp036

####----------------- Section 4: Full reconstruction -- Start -----------------
for ii in range(num_chunk):
    if ii == 0:
        sliceStart = offset + ii*chunk_size
        sliceEnd = offset + (ii+1)*chunk_size
    else:
        sliceStart = offset + ii*(chunk_size-margin_slices)
        sliceEnd = offset + sliceStart + chunk_size
        if sliceEnd > (offset+numRecSlices):
            sliceEnd = offset+numRecSlices

    f = h5py.File(file_name,"r")
    try:
        arr = f["exchange/data"]
    except:
        try:
            arr = f["MyDataDir/TomoShots"]
        except:
            arr = f["entry/data/data"]
    data = arr[0:numProj,sliceStart:sliceEnd,:]
    f.close()
    
    f = h5py.File(flat_name,"r")
    try:
        arr = f["exchange/data_white"]
    except:
        try:
            arr = f["MyDataDir/TomoShots"]
        except:
            arr = f["entry/data/data"]
    white = arr[1:9,sliceStart:sliceEnd,:]
    f.close()
    
    f = h5py.File(dark_name,"r")
    try:
        arr = f["exchange/data_dark"]
    except:
        try:
            arr = f["MyDataDir/TomoShots"]
        except:
            arr = f["entry/data/data"]
    dark = arr[1:9,sliceStart:sliceEnd,:]
    f.close()
    data_size = data.shape
    theta = np.linspace(0,np.pi,num=data_size[0]) 
    print 'data is read'
    
    # replace corrupted images if there is any    
    data[0,:,:] = data[1,:,:]
    
    # remove zingers (pixels with abnormal counts)
    data = tomopy.misc.corr.remove_outlier(data,zinger_level,size=15,axis=0)
    white = tomopy.misc.corr.remove_outlier(white,zinger_level,size=15,axis=0)
    print  'remove outlier is done'
    
    # normalize projection images; for now you need to do below two operations in sequence
    data = tomopy.prep.normalize.normalize(data,white,dark)
    data = tomopy.prep.normalize.normalize_bg(data,air=10)
    print 'normalization is done'
    
    # remove stripes in sinograms
    tomopy.set_debug('True')
    data = tomopy.prep.stripe.remove_stripe_fw(data,level=9,wname='sym16',sigma=2,pad=True)
    tomopy.set_debug('False')
    print 'stripe removal is done'
    
    # phase retrieval
    tomopy.set_debug('True')
    data = tomopy.prep.phase.retrieve_phase(data,pixel_size=pxl,dist=z,energy=eng,alpha=rat,pad=True)
    tomopy.set_debug('False')
    print 'phase retrieval is done'

    # tomo reconstruction
    data = minus_log(data)
    data_recon = recon(data,theta,center=center,algorithm='gridrec')
    print 'reconstruction is done'
    
    # save reconstructions
    tomopy.io.writer.write_tiff_stack(data_recon[np.int(margin_slices/2):(sliceEnd-sliceStart-np.int(margin_slices/2)),:,:], 
                                                 axis = 0,
                                                 fname = output_file, 
                                                 start = sliceStart+np.int(margin_slices/2),
                                                 overwrite = True)
    print 'reconstruction is saved'                                                 
    print time.asctime()                                                  
########----------------- Section 4: Full reconstruction -- End -----------------
