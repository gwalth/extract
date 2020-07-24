#!/usr/bin/env python

#import glob,os,sys
import argparse, os, sys

import astropy.io.fits as pyfits
import numpy as np
from scipy import interpolate
#from scipy.stats.mstats import mquantiles
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties

from math_utils import div_zero
from functions import continuum_fit
from simpledb import simpledb

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike

plot_slit = 0
write = 0

###################################################################################
# usage:
#    #optical_fluxcal.py -dat fgd108.dat --std-src hst
#    #optical_fluxcal.py -dat fgd108.dat --std-src hst -std gd108_LDSS3_1dspec.fits
#    #optical_fluxcal.py -dat fgd108.dat --std-src hst -std GD108_n2/gd108_LDSS3_1dspec.fits -std-og GD108_OG_n2/gd108_LDSS3_OG_1dspec.fits
#    optical_fluxcal.py -dat fgd108.dat --std-src hst -std GD108_n2/gd108_LDSS3_1dspec.fits -std-og GD108_OG_n2/gd108_LDSS3_OG_1dspec.fits -dir macs1314/macs1314 -smf macs1314/macs1314.SMF
#
#
#    cd /Users/gwalth/data/LDSS3/REDUCTIONS/kelson
#    optical_fluxcal.py -dat fltt377.dat --std-src ctio -std ut180712/LTT377/LTT377_center_1dspec.fits  -dir HLS/macs1314/macs1314 -smf HLS/macs1314/macs1314.SMF --external-response ldss3_center_reponse_crap.fits
###################################################################################

parser = argparse.ArgumentParser(description='flux calibrate 1d spectra')

parser.add_argument('-sci', metavar='file', type=str, nargs='?',
                    default=None,help='raw science FITS file')
parser.add_argument('-std', metavar='file', type=str, nargs='?',
                    default=None,help='raw standard FITS file')
parser.add_argument('-std-og', metavar='file', type=str, nargs='?',
                    default=None,help='raw standard OG590 filter FITS file')
parser.add_argument('-dat', metavar='file', type=str, nargs='?', required=True,
                    default=None,help='standard file photometry/flux file')
parser.add_argument('-fit', choices=['spline','legendre'], default="spline",
                    help='fit to flux standard') 
parser.add_argument('-ord', metavar='order of fit', type=int, nargs='?',
                    default=5, help='order of fit')
parser.add_argument('--std-src', choices=['ctio','hst','snf','oke'], default="ctio",
                    help='standard template set') 

parser.add_argument('--external-response', metavar='file', type=str, nargs='?',
                    default=None,help='Use external reponse function to flux calibrate')
parser.add_argument('-smf', metavar='file', type=str, nargs='?',
                    default=None,help='SMF file')
parser.add_argument('-dir', metavar='file', type=str, nargs='?',
                    default=None,help='directory of 2d spectra')


args = parser.parse_args()

std     = args.std
std_og  = args.std_og
sci     = args.sci
dat     = args.dat
std_src = args.std_src
fit     = args.fit
order   = args.ord

external_response = args.external_response


smf  = args.smf
fdir = args.dir


ph1 = os.getenv("PYTHONHOME1")
if not ph1: raise "PYTHONHOME1 environmental variable not set!"

std_dir = ph1 + "/datafiles/" + std_src 
if std_src != "snf": std_dir += "stan"

dat_file = std_dir + '/' + dat

if not os.path.exists(dat_file):
    print
    print "Standard %s not found!!!" % (dat)
    print
    print "%s available standards:" % (std_src)
    print os.listdir(std_dir)
    
    sys.exit()

data = np.loadtxt(dat_file)
print data

wav_std = data[:,0]   # column 0
flx_std = data[:,1]*10**-16   # column 1  ergs/cm^2/s/A    F_lambda
noi_std = data[:,2]*10**-16   # column 2  ergs/cm^2/s/A    F_lambda (noise)
#flx_std = data[:,1]   # column 1 [10^-16] ergs/cm^2/s/A    F_lambda
#noi_std = data[:,2]   # column 2 [10^-16] ergs/cm^2/s/A    F_lambda (noise)
#fname = f[1:].replace(".dat","")

ext = 0
# read in frame to be used as a standard
pf = pyfits.open(std)
print len(pf)
print pf.info
head = pf[ext].header
naxis1   = head["naxis1"]
naxis2   = head["naxis2"]
crval1   = head["crval1"]
cdelt1   = head["cdelt1"]
crpix1   = head["crpix1"]
dcflag   = head["dc-flag"]
exptime  = head["exptime"]
airmass  = head["airmass"]

stdflt_exptime = head["fltexp"] # std flat exposure

wav_obj = (np.arange(naxis1)+1-crpix1)*cdelt1+crval1
if dcflag: wav_obj = np.power(10,wav_obj)

data = pf[ext].data

print data.shape

flx_obj = data[0,::]
noi_obj = data[1,::]
flt_obj = data[3,::]

print flx_obj.shape
print noi_obj.shape

#print stdflt_exptime

nflx_obj = flx_obj/exptime
nnoi_obj = noi_obj/exptime
#nflt_obj = flt_obj/(np.median(flt_obj)*stdflt_exptime)
nflt_obj = flt_obj/(stdflt_exptime)
gain1 = np.median(nflt_obj)

print gain1


# ctio extinction
# /iraf/iraf/noao/lib/onedstds
cef = "ctioextinct.dat"
data = np.loadtxt(ph1+"/datafiles/"+cef)
wav_ext = data[:,0]   # column 0
mag_ext = data[:,1]   # column 1  magnitudes

print wav_ext
# interpolation of extinction
Ie = interpolate.interp1d(wav_ext,mag_ext,bounds_error=False,fill_value=0)
Ie_ext = extrap1d(Ie)

if std_og:
    pf = pyfits.open(std_og)
    print len(pf)
    print pf.info
    head = pf[ext].header
    naxis1   = head["naxis1"]
    naxis2   = head["naxis2"]
    crval1   = head["crval1"]
    cdelt1   = head["cdelt1"]
    crpix1   = head["crpix1"]
    dcflag   = head["dc-flag"]
    airmass  = head["airmass"]

    exptime_og  = head["exptime"]

    stdflt_exptime_og = head["fltexp"] # std flat exposure

    wav_obj_og = (np.arange(naxis1)+1-crpix1)*cdelt1+crval1
    if dcflag: wav_obj_og = np.power(10,wav_obj_og)

    data = pf[ext].data

    print data.shape

    flx_obj_og = data[0,::]
    noi_obj_og = data[1,::]
    flt_obj_og = data[3,::]

    print flx_obj_og.shape
    print noi_obj_og.shape


    nflx_obj_og = flx_obj_og/exptime_og
    nnoi_obj_og = noi_obj_og/exptime_og
    #nflt_obj = flt_obj/(np.median(flt_obj)*stdflt_exptime)
    nflt_obj_og = flt_obj_og/(stdflt_exptime_og)
    gain2 = np.median(nflt_obj_og)

    print gain2



#print wav_std
#print wav_obj
# interpolation of the flux standard 
Is = interpolate.interp1d(wav_std,flx_std,bounds_error=False,fill_value=0)


if wav_std[0] < wav_ext[0]: w1 = wav_ext[0]
else: w1 = wav_std[0]
if wav_std[-1] > wav_ext[-1]: w2 = wav_ext[-1]
else: w2 = wav_std[-1]

filt = np.ones((wav_obj.shape))
if wav_obj[0] < w1:
   filt = wav_obj > w1
if wav_obj[-1] > w2:
   filt = wav_obj < w2
        


# remove points beyond interpolation region
wav_obj  = np.compress(filt,wav_obj,0)
flx_obj  = np.compress(filt,flx_obj,0)
nflx_obj = np.compress(filt,nflx_obj,0)
nnoi_obj = np.compress(filt,nnoi_obj,0)
nflt_obj = np.compress(filt,nflt_obj,0)
print len(wav_obj)


UV       = 1.0*np.less(wav_obj,4200)
Hepsilon = 1.0*np.greater(wav_obj,4065)*np.less(wav_obj,4150)
Hgamma   = 1.0*np.greater(wav_obj,4300)*np.less(wav_obj,4400)
Hbeta    = 1.0*np.greater(wav_obj,4765)*np.less(wav_obj,4920)
Halpha   = 1.0*np.greater(wav_obj,6520)*np.less(wav_obj,6630)
Bband    = 1.0*np.greater(wav_obj,6850)*np.less(wav_obj,6950)
unknown  = 1.0*np.greater(wav_obj,7100)*np.less(wav_obj,7350)
Aband    = 1.0*np.greater(wav_obj,7580)*np.less(wav_obj,7700)
#
Balmer = Hepsilon + Hgamma + Hbeta + Halpha
Atmos = Aband + Bband
#
#bad = Atmos+Balmer+UV
bad = Atmos+Balmer+unknown
good = np.logical_not(bad)

# atmospheric extinction correction for std
extinction = Ie_ext(wav_obj)
e_nflx_obj = nflx_obj*10**(0.4*extinction*airmass)
e_nnoi_obj = nnoi_obj*10**(0.4*extinction*airmass)


# flux / spectrophotometric standard

# response of the standard spectrum
print wav_obj
print len(wav_obj)
R = np.array([div_zero(e_nflx_obj[j],Is(wo)) for j,wo in enumerate(wav_obj)])
N = np.array([div_zero(e_nnoi_obj[j],Is(wo)) for j,wo in enumerate(wav_obj)])
#R = np.array([div_zero(nflx_obj[j],Is(wo)) for j,wo in enumerate(wav_obj)])
#N = np.array([div_zero(nnoi_obj[j],Is(wo)) for j,wo in enumerate(wav_obj)])
print len(R)
print len(N)


weight = div_zero(good.astype(np.float32),N)
#print len(weight)

#print R
#print wav_obj
#print weight
print
print
print
print

def spline_fit1(wav,spec,weight,N2=20):
    #print 
    ssize = (np.max(wav) - np.min(wav))/N2
    #print ssize
    wav_knot = np.arange(np.min(wav),np.max(wav),ssize)[1:]
    print spec
    #print wav
    #print wav_knot
    #print weight
    print len(spec)
    #print len(wav)
    #print len(wav_knot)
    #print len(weight)

    # scipy 0.14
    #C = interpolate.LSQUnivariateSpline(wav,spec,wav_knot,w=weight)
    # scipy 1.5
    C = interpolate.LSQUnivariateSpline(wav,spec,wav_knot,w=weight,ext=1)
    # ext=1 extrapolation mode return zeros

    spec_fn = C.__call__(wav)
    spec_knot = C.__call__(wav_knot)
    #spec_fn = C(wav)
    #spec_knot = C(wav_knot)
    print C.get_knots()
    print "Rchi2 =",C.get_residual()

    print spec_knot
    print
    print
    print
    return spec_fn,spec_knot,wav_knot

def spline_fit2(wav,spec,weight):
    #print 
    #C = interpolate.UnivariateSpline(wav,spec,w=weight,k=3,s=4)
    C = interpolate.UnivariateSpline(wav,spec,w=weight,k=3,s=4,ext=1)
    print spec
    print wav
    print weight
    spec_fn = C.__call__(wav)
    wav_knot = C.get_knots()
    spec_knot = C.__call__(wav_knot)
    print wav_knot
    #spec_fn = C(wav)
    #spec_knot = C(wav_knot)
    print C.get_knots()
    print "Rchi2 =",C.get_residual()

    print spec_knot
    print
    print
    print
    return spec_fn,spec_knot,wav_knot


##################################
## spline fit to response function
##################################
#R_fn,R_dict = continuum_fit(R,wav_obj,weight,fit="spline",full_output=1,N2=15)
#wav_knot = R_dict["wav_knot"]
#R_knot   = R_dict["spec_knot"]

R_fn, R_knot, R_wav_knot = spline_fit1(wav_obj, R, weight, N2=16)
#R_fn, R_knot, R_wav_knot = spline_fit1(wav_obj, R, weight, N2=30)
#R_fn, R_knot, R_wav_knot = spline_fit2(wav_obj, R, weight)

print "spline fit to response function"
print np.max(R_knot)
print np.max(R_fn)
print
print



######################
## spline fit to flat
######################
print "spline fit to flat"
#print nflt_obj
#print wav_obj
#print weight
#print
#F_fn,F_dict = continuum_fit(nflt_obj,wav_obj,weight,fit="spline",full_output=1,N2=15)
#wav_knot = F_dict["wav_knot"]
#F_knot   = F_dict["spec_knot"]


F_fn, F_knot, F_wav_knot = spline_fit1(wav_obj, nflt_obj, weight, N2=16)
#F_fn, F_knot, F_wav_knot = spline_fit2(wav_obj, nflt_obj, weight)
rawflat_max = np.max(F_knot)
print rawflat_max
print np.max(F_fn)
print
print



# slit to slit correction 
nR_fn = div_zero(R_fn,F_fn/rawflat_max)

w0 = wav_obj[0]
w1 = wav_obj[-1]

fig = plt.figure()
p1 = fig.add_subplot(111)
p1.plot(wav_obj,nR_fn,c='k')
#p1.set_title("%s Spectra" % fname)
p1.set_title("Flat Corrected Response")
p1.set_xlabel("Wavelength ($\AA$)")
p1.set_xlim(w0,w1)
plt.show()




#####################
# check standard fits
#####################
font = FontProperties()
font.set_size(10)


fig = plt.figure()

w0 = wav_obj[0]
w1 = wav_obj[-1]

y0 = np.min(flx_std)
y1 = np.max(flx_std)

p1 = fig.add_subplot(111)
#p1 = fig.add_subplot(111)
p1.plot(wav_obj,nflx_obj,c='g',label="Normalized Sci")
if std_og:
    p1.plot(wav_obj_og,nflx_obj_og,c='g',alpha=0.5,label="Normalized Sci OG")
#p1.set_title("%s Spectra" % fname)
#p1.set_title("Raw spectra")
p1.set_xlabel("Wavelength ($\AA$)")
p1.set_ylabel("Flux (counts/s)")
p1.set_title(std)
p1.set_xlim(w0,w1)

leg = p1.legend(loc=1,numpoints=1,prop=font)
leg.draw_frame(False)

plt.show()



fig = plt.figure(figsize=(10,10))

w0 = wav_obj[0]
w1 = wav_obj[-1]

y0 = np.min(flx_std)
y1 = np.max(flx_std)

p1 = fig.add_subplot(221)
#p1 = fig.add_subplot(111)
p1.plot(wav_std,flx_std,c='k',label="Reference")
p1.plot(wav_obj,flx_obj,c='r',label="Raw Sci")
p1.plot(wav_obj,nflx_obj,c='g',label="Normalized Sci")
p1.plot(wav_obj,nflt_obj,c='c',label="Raw Flat")
if std_og:
    p1.plot(wav_obj_og,flx_obj_og,c='r',alpha=0.5,label="Raw Sci OG")
    p1.plot(wav_obj_og,nflx_obj_og,c='g',alpha=0.5,label="Normalized Sci OG")
    p1.plot(wav_obj_og,nflt_obj_og,c='c',alpha=0.5,label="Raw Flat OG")
#p1.set_title("%s Spectra" % fname)
p1.set_title("Raw spectra")
p1.set_xlabel("Wavelength ($\AA$)")
p1.set_xlim(w0,w1)

leg = p1.legend(loc=1,numpoints=1,prop=font)
leg.draw_frame(False)

p2 = fig.add_subplot(222)
p2.plot(wav_obj,R,c='k',label="Response")
p2.plot(wav_obj,N,c='r')
p2.plot(wav_obj,R_fn,c='g',label="Fit to Response")
if fit == "spline":
    p2.scatter(R_wav_knot,R_knot,c='b')
#p2.set_title("%s Response" % fname)
p2.set_title("Response function")
p2.set_xlabel("Wavelength ($\AA$)")
p2.set_xlim(w0,w1)

#leg = p2.legend(loc=1,numpoints=1,prop=font)
#leg.draw_frame(False)

y1 = 1.2*np.max(F_knot)

p3 = fig.add_subplot(223)
#p3.plot(wav_obj,nflt_obj,c='k',label="Raw Flat")
p3.plot(wav_obj,F_fn,c='g',label="Fit to Raw Flat")
if fit == "spline":
    p3.scatter(F_wav_knot,F_knot,c='b')
#p3.set_title("%s Flat" % fname)
p3.set_title("Fit to raw spectral flat")
p3.set_xlabel("Wavelength ($\AA$)")
p3.set_ylabel("Flux (counts)")
p3.set_xlim(w0,w1)
p3.set_ylim(0,y1)

#leg = p3.legend(loc=1,numpoints=1,prop=font)
#leg.draw_frame(False)

corrected = e_nflx_obj/R_fn
y1 = np.max(1.2*corrected)

p4 = fig.add_subplot(224)
p4.plot(wav_obj,corrected,c='r',label="Corrected")
p4.plot(wav_std,flx_std,c='k',label="Reference")
p4.set_xlabel("Wavelength ($\AA$)")
p4.set_ylabel("F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2} \AA^{-1}$")
#p4.set_title(fname)
p4.set_title("Flux calibrated spectra")
p4.set_xlim(w0,w1)
#p4.set_ylim(y0,y1)
p4.set_ylim(0,y1)

#leg = p4.legend(loc=1,numpoints=1,prop=font)
#leg.draw_frame(False)

if std_og:
##############################
# derive OG filter emperically
##############################

    fig = plt.figure()

    ymin    = np.min(nflt_obj[0:100])
    ymin_og = np.min(nflt_obj_og[0:100])
    print ymin
    
    y_num = np.where(np.less_equal(nflt_obj_og - ymin_og,0),0,nflt_obj_og - ymin_og)
    y_den = np.where(np.less_equal(nflt_obj - ymin,0),0,nflt_obj - ymin)

    og_filter = y_num / y_den

    p1 = fig.add_subplot(121)
    p1.plot(wav_obj,nflt_obj,c='c',label="Raw Flat")
    p1.plot(wav_obj_og,nflt_obj_og,c='c',alpha=0.5,label="Raw Flat OG")
    p1.set_title("Raw spectra")
    p1.set_xlabel("Wavelength ($\AA$)")
    p1.set_xlim(w0,w1)

    leg = p1.legend(loc=2,numpoints=1,prop=font)
    leg.draw_frame(False)

    p2 = fig.add_subplot(122)
    p2.plot(wav_obj,og_filter,c='k',label="OG filter curve")
    #p2.plot(wav_obj,nflt_obj_og/nflt_obj,c='c',label="OG filter curve")
    p2.set_xlabel("Wavelength ($\AA$)")
    p2.set_title("Empirical OG590 filter")
    p2.set_xlim(w0,w1)

    #plt.show()
    #sys.exit()



##############
# read in mask
##############


db = simpledb(smf,fdir)
#print db
objD = db.dict["objects"]
rows = db.row_search()
print rows
for d in objD: print d,objD[d]

N = len(rows)


smf_x = [objD[rows[fr]]["setup"]["smf_x"] for fr in xrange(N) if objD[rows[fr]]["setup"]["smf_type"] == "SLIT"]
smf_y = [objD[rows[fr]]["setup"]["smf_y"] for fr in xrange(N) if objD[rows[fr]]["setup"]["smf_type"] == "SLIT"]

#NN = len(smf_x)
#colors = [cmap(i) for i in np.linspace(0, 0.9, NN)]

cmap = cm.gist_ncar
#cmap = cm.jet


# y = m*x + c
# 0 = m*vmin + c --> c = -m*vmin
# 1 = m*vmax + c --> m = (1-c)/vmax --> m = 1-(-m*vmin)/vmax -->   m = (1 + m*vmin)/vmax --> vmax*m = 1 + vmin*m -->  m(vmax-vmin)=1
#                    m = 1/(vmax-vmin)
#colors = [cmap(x) for x in smf_x]

vmin = np.min(smf_x)
vmax = np.max(smf_x)
colors = [cmap( x/(vmax-vmin) - vmin/(vmax-vmin) ) for x in smf_x]
##colors = [x/(vmax-vmin) - vmin/(vmax-vmin)  for x in smf_x]


vmin = np.min(smf_y)
vmax = np.max(smf_y)
colors = [cmap( y/(vmax-vmin) - vmin/(vmax-vmin) ) for y in smf_y]
##colors = [y/(vmax-vmin) - vmin/(vmax-vmin)  for y in smf_y]
#print
#print
#print
#print smf_x
#print
#print smf_y
#print
#print colors
#print
#print
#print









#fig, ax = plt.subplots(figsize=(6, 1))
#fig.subplots_adjust(bottom=0.5)

#norm = mpl.colors.Normalize(vmin=-100, vmax=100)

#cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
##cb1 = mpl.colorbar.ColorbarBase(p1, cmap=cmap, norm=norm, orientation='horizontal')
##cb1 = mpl.colorbar.ColorbarBase(p1, cmap=cmap, norm=norm, orientation='horizontal')
#cb1.set_label('Mask position')
##plt.show()


if std_og:
##########################################
    # splice OG and non-OG filter
    match_min = 6600
    match_max = 6850

    splice_cen = 6725

    # non-OG
    match_min_ind = np.argmin(abs(wav_obj - match_min))
    match_max_ind = np.argmin(abs(wav_obj - match_max))
    #print match_min_ind
    #print match_max_ind
    #print wav_obj[match_min_ind]
    #print wav_obj[match_max_ind]
    #print (wav_obj[match_min_ind] + wav_obj[match_max_ind])/2
    med_bp = np.median(nflx_obj[match_min_ind:match_max_ind])
    #print med_bp
    #print
    
    # OG
    match_min_ind = np.argmin(abs(wav_obj_og - match_min))
    match_max_ind = np.argmin(abs(wav_obj_og - match_max))
    #print match_min_ind
    #print match_max_ind
    #print wav_obj_og[match_min_ind]
    #print wav_obj_og[match_max_ind]
    
    #print (wav_obj_og[match_min_ind] + wav_obj_og[match_max_ind])/2


    med_bp_og = np.median(nflx_obj_og[match_min_ind:match_max_ind])
    #print med_bp_og


    offset = med_bp - med_bp_og
    #print offset


    # blue non-og + red og
    #   < 6720        > 6720 + offset
    cen_ind = np.argmin(abs(wav_obj - splice_cen))
    #print cen_ind
    
    #print nflx_obj.shape
    #print nflx_obj_og.shape


    nflx_obj_com = nflx_obj.copy()
    #print nflx_obj_com.shape


    # composite spectrum, need to test with lots of spectra!!!
    #nflx_obj_com[cen_ind:-1] = nflx_obj_og[cen_ind:-1]/og_filter[cen_ind:-1] - offset
    nflx_obj_com[cen_ind:-1] = nflx_obj_og[cen_ind:-1] + offset
    #nflx_obj_com[cen_ind:-1] = nflx_obj_og[cen_ind:-1]



    fig = plt.figure(figsize=(10,10))
    
    w0 = wav_obj[0]
    w1 = wav_obj[-1]
    
    y0 = np.min(flx_std)
    y1 = np.max(flx_std)
    
    p1 = fig.add_subplot(111)
    #p1 = fig.add_subplot(111)
    p1.plot(wav_obj,nflx_obj,c='g',label="Normalized Sci")
    p1.plot(wav_obj_og,nflx_obj_og,c='g',alpha=0.5,label="Normalized Sci OG")
    p1.plot(wav_obj,nflx_obj_com,c='r',label="Composite Sci")
    #p1.set_title("%s Spectra" % fname)
    #p1.set_title("Raw spectra")
    p1.set_xlabel("Wavelength ($\AA$)")
    p1.set_xlim(w0,w1)
    
    leg = p1.legend(loc=1,numpoints=1,prop=font)
    leg.draw_frame(False)

    plt.show()

    #wav_obj,nflx_obj
    #wav_obj_og,nflx_obj_og
    
    e_nflx_obj_com = nflx_obj_com*10**(0.4*extinction*airmass)
    
    # flux / spectrophotometric standard
    
    # response of the standard spectrum
    R_com = np.array([div_zero(e_nflx_obj_com[j],Is(wo)) for j,wo in enumerate(wav_obj)])
    
    
    
    R_min_ind = np.argmin(abs(wav_obj_og - 7000))
    R_max_ind = np.argmin(abs(wav_obj_og - 8000))
    R_max = np.max(R_com[R_min_ind:R_max_ind])
    
    R_com /= R_max # normalize
    
    R_fn, R_knot, R_wav_knot = spline_fit1(wav_obj, R_com, weight, N2=17)
    
    
    fig = plt.figure(figsize=(10,10))
    p1 = fig.add_subplot(111)
    #p1.plot(wav_obj,R,c='k',label="Response")
    p1.plot(wav_obj,R_com,c='g',label="Response Composite")
    p1.plot(wav_obj,R_fn,c='b',label="Fit to Response")
    if fit == "spline":
        p1.scatter(R_wav_knot,R_knot,c='b')
    #p1.set_title("%s Response" % fname)
    p1.set_title("Response function")
    p1.set_xlim(w0,w1)
    
    leg = p1.legend(loc=1,numpoints=1,prop=font)
    leg.draw_frame(False)
    
    plt.show()
    
    hdu = pyfits.PrimaryHDU([R_fn,wav_obj],header=header)
    hdu.writeto("ldss3_center_reponse_crap.fits",clobber=True)
    
    
    sys.exit()
    if write:
        plt.savefig("%s_FluxCalib.ps" % (psname))
    else:
    
        plt.show()


if external_response:
    pf = pyfits.open(external_response)
    print len(pf)
    print pf.info
    head = pf[ext].header
    naxis1   = head["naxis1"]
    naxis2   = head["naxis2"]
    crval1   = head["crval1"]
    cdelt1   = head["cdelt1"]
    crpix1   = head["crpix1"]
    dcflag   = head["dc-flag"]
    exptime  = head["exptime"]
    airmass  = head["airmass"]

    data = pf[ext].data

    print data.shape

    norm_R = data[0,::]
    wav_norm_R = data[1,::]

    fig = plt.figure()

    p1 = fig.add_subplot(111)
    p1.plot(wav_norm_R,norm_R,c='k',label="External 2nd Order Corrected Response")
    p1.plot(wav_obj,R_fn/np.max(R_fn),c='g',label="Response From Star")
    p1.set_title("Response function")
    p1.set_xlim(w0,w1)
    p1.set_ylim(-0.1,1.25)

    leg = p1.legend(loc=1,numpoints=1,prop=font)
    leg.draw_frame(False)

    I_norm_R = interpolate.interp1d(wav_norm_R,norm_R,bounds_error=False,fill_value=0)
    
    plt.show()

#I_R = interpolate.interp1d(wav_obj,R_fn,bounds_error=False,fill_value="extrapolate")
I_R = interpolate.interp1d(wav_obj,R_fn,bounds_error=False,fill_value=0)
#I_R = interpolate.interp1d(wav_obj,nR_fn,bounds_error=False,fill_value=0)
print wav_obj
print wav_obj[0]
print wav_obj[-1]
#sys.exit()


#fig = plt.figure()
#p1 = fig.add_subplot(111)


i = 0
for fr in xrange(N):
    setup_obj = objD[rows[fr]]["setup"]
    fspec = setup_obj["spec1d"]
    obj = setup_obj["object"]
    smftype = setup_obj["smf_type"]
    print "row = %i" % fr
    print "obj = %s" % obj
    print fspec
    #print

    #if setup_obj["smf_type"] == "SLIT":
    #or setup_obj["smf_type"] == "HOLE":
    if 1:

        try: 
            pf = pyfits.open(fspec,ignore_missing_end=1)
        except:
            print "something went wrong reading!!!"
            print
            print
            print
            continue

        spec1d = pf[0].data
        head = pf[0].header
      
        crpix1  = head['CRPIX1']     # starting pixel
        crval1  = head['CRVAL1']     # starting wavelength
        cd1_1   = head['CD1_1']      # dispersion
        dc_flag = head['DC-FLAG']    # Log-linear flag

        naxis1  = head["naxis1"]
        exptime = head["exptime"]
        airmass = head["airmass"]

        #sciflt_exptime = head["fltexp"] # sci flat exposure
        print spec1d.shape
        
        ysize,xsize = spec1d.shape
        x = np.arange(xsize)+1
        wav = (x-crpix1)*cd1_1+crval1  # compute wavelength
        
        if dc_flag: wav = np.power(10,wav)

        sci = spec1d[0,::]
        noi = spec1d[1,::]
        flt = spec1d[3,::]

     





        print wav.shape
        #print sci.shape
        #print noi.shape

        print np.max([wav[0],wav_std[0],wav_ext[0]])
        print np.min([wav[-1],wav_std[-1],wav_ext[-1]])


        #if wav_std[0] < wav_ext[0]: w1 = wav_ext[0]
        #else: w1 = wav_std[0]
        #if wav_std[-1] > wav_ext[-1]: w2 = wav_ext[-1]
        #else: w2 = wav_std[-1]

        #print wav_std[0],wav_ext[0],wav[0]
        #print wav_std[-1],wav_ext[-1],wav[-1]
        #print w1,w2

        #filt = np.ones((wav.shape))
        #if wav[0] < w1:
        #   filt = wav > w1
        #if wav[-1] > w2:
        #   filt = wav < w2
        
        #print w1,w2
        
        #print "sci",wav
        #print "ext",wav_ext
        #print "std",wav_std
        ## remove points beyond interpolation region
        #wav  = np.compress(filt,wav,0)
        #sci  = np.compress(filt,sci,0)
        #noi  = np.compress(filt,noi,0)
        #flt  = np.compress(filt,flt,0)
        #print len(wav)
        #print wav.shape

        nsci = sci/exptime
        nnoi = noi/exptime
        #nflt = flt/(sciflt_exptime)
        nflt = flt/np.max(flt)

        print wav
        extinction2 = Ie_ext(wav)
        R2  = I_R(wav)  # interpolated response
        #flx = Is(wav)   # interpolate flux from standard
        if external_response:
            R2 = R2*I_norm_R(wav) # interpolated normalized response




        #print extinction2.shape
        #print R.shape



        e_nsci = nsci*10**(0.4*extinction2*airmass)
        e_nnoi = nnoi*10**(0.4*extinction2*airmass)
        final_sci = e_nsci/R2  # original 
        final_noi = e_nnoi/R2
        #p1.plot(wav,flt,color=colors[i]) # plot raw flats
        final_sci_nflt = final_sci/nflt
        final_noi_nflt = final_noi/nflt


        if plot_slit:
            fig2 = plt.figure()
            p1 = fig2.add_subplot(121)
            p1.plot(wav,nsci,c='k',label="science")
            p1.plot(wav,nnoi,c='r',label="noise")
            p1.plot(wav,nflt,c='b',label="flat")
            p1.set_yscale("log")
            p1.set_ylabel("Flux (counts/s)")
            p1.set_xlabel("Wavelength ($\AA$)")
            leg1 = p1.legend(loc=2,numpoints=1,prop=font)
    
            p2 = fig2.add_subplot(122)
            p2.plot(wav,final_sci,c='k',label="science")
            p2.plot(wav,final_sci_nflt,c='k',label="science / normalized flat",alpha=0.5)
            p2.plot(wav,final_noi,c='r',label="noise")
            p2.plot(wav,final_noi_nflt,c='r',label="noise / normalized flat",alpha=0.5)
            p2.set_yscale("log")
            p2.set_xlabel("Wavelength ($\AA$)")
            p2.set_ylabel("F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2} \AA^{-1}$")
            leg2 = p2.legend(loc=2,numpoints=1,prop=font)
    
            fig2.suptitle("%s   row = %i   id = %s   N = %i" % (smftype,fr,obj,N))
    
            plt.show()

        #try: 
        #new_data = [final_sci,final_noi] # for testing

        new_head = head.copy()
        new_head["ARRAY7"] = ("FLUXED SPECTRUM","units of erg/s/cm^2/Ang")
        new_head["ARRAY8"] = ("FLUXED NOISE","units of erg/s/cm^2/Ang")

        #new_data = spec1d + [final_sci,final_noi]
        new_data = np.concatenate([spec1d,final_sci[np.newaxis,::],final_noi[np.newaxis,::]],0)
        print new_data.shape
        hdu = pyfits.PrimaryHDU(new_data,header=new_head)
        #hdu = pyfits.PrimaryHDU([final_sci,final_noi,wav],header=head)
        #hdu = pyfits.PrimaryHDU([final_sci_nflt,final_noi_nflt,wav],header=head)
        nfspec = fspec.replace(".fits","_flux.fits")
        hdu.writeto(nfspec,clobber=True)
        #except:
        #    print "something went wrong writing!!!"
        #    print
        #    print
        #    print
        #    continue

        # method:
        #   flux_sci/counts_sci = flux_star/counts_star
        #   flux_sci = flux_star/counts_star*counts_sci
        #   R = counts_star/flux_star
        #   flux_sci = counts_sci/R

        
        #p2.plot(wav_std,flx_std,c='k',label="Reference")

        i+= 1

        print
        print
        print

#cb = mpl.colorbar.ColorbarBase(p1, cmap=cmap)
#cb = fig.colorbar.ColorbarBase(p1, cmap=cmap)
#cb = mpl.colorbar(p1, cmap)


#cb = fig.colorbar(color)
#cb = fig.colorbar(px)
#cb = fig.colorbar(p1,boundaries=arange(-4000,4000+50,50))
#p1.set_title("Raw spectral flats (color-coded in mask pos in X)")

#p1.set_title("Raw spectral flats (color-coded in mask pos in Y)")
#p1.set_xlabel("Wavelength ($\AA$)")
#p1.set_ylabel("Counts")
#plt.show()

