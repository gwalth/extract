#!/usr/bin/env python

######################################
# From Dan's code
#from basics import *
#from VTKHelperFunctions import *
#from VTKHelperFunctions import VTKDivide, VTKMultiply, VTKGreaterEqual,VTKGauss
######################################
from basics import *
from FITS import *
from Multipack import *
from VTKHelperFunctions import *

import argparse,sys
import inspect
import cPickle

import pyfits
import numpy as np
from scipy import signal
from scipy.optimize import leastsq
import numpy.linalg as linalg

import matplotlib.pyplot as plt

from simpledb import simpledb

# https://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy
#def consecutive(data, stepsize=1):
#    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

#consecutive(a)

def between(x,x0,x1):
    return np.greater(x,x0)*np.less(x,x1)

def legendre(x,nl):
    x = x.astype("Float64")
    l = (nl+1,)+x.shape
    P = np.zeros(l,dtype=np.float64)
    P[0] = np.ones(x.shape,dtype=np.float64)
    if nl > 0: P[1] = x
    if nl > 1:
       for n in range(1,nl):
           P[n+1] = ((2*n+1)*x*P[n] - n*P[n-1])/(n+1)
    return P


def param(profile,x,xcen,verb=0):
    x0,x1 = np.clip([xcen-3,xcen+4],0,x[-1])
    wt = np.greater(profile,0)
    wt2 = np.asarray([wt,wt])
    kz, = np.nonzero(wt)

    #if kz:
    # this part is not functional
    wstatus = 1
    if 0:
        wstatus = 1
        j = kz[np.argmin(np.abs(kz-xcen))]
        print j
        #wt = 1.0*VTKSeedConnectivity(wt2,x=[j],y=[0])[0]
        print wt2
        wt = consecutive(wt2)
        print wt
        sys.exit()

        #wt = VTKDilate(wt,5,numret=1)
        #qdump("wt.fits",1.0*wt)
        #qdump("swt.fits",1.0*swt)

        # counting the number of weights, should not be a small number!
        if sum(wt) > 2:
            if verb: print "profile",profile
            if verb: print "xcen",xcen
            if verb: print "wt",wt
            xcen = divz(sum(x*wt),sum(wt))
            xcen = int(xcen+0.5)
            if verb: print "xcen",xcen
            x0,x1 = np.clip([xcen-3,xcen+4],0,x[-1])
        else: wstatus = 0

    ### wt fails kz is a zero array!!!
    #else: wstatus = 0

    if x1-x0 < 7:
       if x0 == 0: x1 = 7
       if x1 == x[-1]: x0 = x1-7
    j = np.argmax(profile[x0:x1]) + x0
    xcen = x[j]
    amp = profile[j]
    if amp > 0: amp = np.log10(amp)
    else: amp = 2.0
    sig = np.log10(1.5)
    #inside = between(x,x0,x1)
    inside = np.greater(x,x0)*np.less(x,x1)
    outside = np.logical_not(inside)
    const = np.median(np.compress(outside,profile))
    if verb: print xcen,amp,sig,const
    return [xcen,amp,sig,const],wstatus

def gauss(param,data,x,verb=0,fitting=1):
    xcen,amp,sig,const = param
    amp = np.power(10,amp)
    sig = np.power(10,sig)
    if verb: print xcen,amp,sig,const
    gfunc = amp*np.exp(-0.5*((x-xcen)/sig)**2.) + const

    if fitting: return data - gfunc
    else: return gfunc

def GetTrace(D,params,edge,trace_order=2,sw=10):
    print D["data"]
    #(spec2d,noise2d) = D["data"]
    spec2d = D["data"]["spec2d"]
    noise2d = D["data"]["noise2d"]
    N = spec2d.shape
    y, x = np.indices(N,dtype=np.float32) 

    #if "LDSS3" in instrument: ReallyBad = less(y,5)
    #else: ReallyBad = logical_or(less(y,2),greater(y,y.shape[0]-3))
    EdgeMask = np.logical_or(np.less(y,edge),np.greater(y,y.shape[0]-(1+edge)))
    y = y - params[0]
    x0, x1 = 0.0, 1.0*x.shape[1]
    xc, xw = (x1+x0)/2.0, (x1-x0)/2.0

    w = divz(spec2d*np.greater_equal(spec2d,0),noise2d)
    w = w*np.greater_equal(noise2d,0)
    w = w/(1+divz(spec2d-signal.medfilt(spec2d,kernel_size=7),noise2d)**2)

    # VTK
    #w = VTKDivide(VTKMultiply(spec2d,VTKGreaterEqual(spec2d,0)),noise2d)
    #w = VTKMultiply(w,VTKGreaterEqual(noise2d,0))
    #w = w/(1+divz(spec2d-VTKMedian(spec2d,7,1,numret=1),noise2d)**2)

    #w = w*greater(spec2d,-3*noise2d)
    #w = w*greater(w,0.0)*logical_not(ReallyBad)
    w = w*np.greater(w,0.0)*np.logical_not(EdgeMask)
    w0 = w*between(y,-sw,+sw)
    cy = divz(np.sum(y.flatten()*w0.flatten()**2),np.sum(w0.flatten()**2))
    w = w*between(y-cy,-sw,+sw)
    b = legendre((x.flatten()-xc)/xw,trace_order)
    t = svdfit(np.transpose(b)*w.flatten()[::,np.newaxis],y.flatten()*w.flatten())
    #sys.exit()
    return t
#    Test the solution by rebinning spectra:
#    p = reshape(dot(t,b),x.shape)
#    r = VTKImageTransform(spec2d,0.0*p,p,reverse=0,numret=1)

def GetProfile(D,edge,plot=0,verb=0):
    #(spec2d,noise2d) = D["data"]
    spec2d = D["data"]["spec2d"]
    noise2d = D["data"]["noise2d"]
    N = spec2d.shape
    y, x = np.indices(N,dtype=np.float32) 
    EdgeMask = np.logical_or(np.less(y,edge),np.greater(y,y.shape[0]-(1+edge))) 
    spec2d_noedge = spec2d*np.logical_not(EdgeMask)                       # w in Dan's code

    xcen = D["pos"]
    obj = D["id"]
    spec1d_noedge = np.sum(spec2d_noedge,axis=0)                          # pz in Dan's code
    spec1d_clip = np.not_equal(spec1d_noedge,0)                           # gz in Dan's code

    #profile1 = np.median(spec2d,axis=1)
    #profile2 = np.median(spec2d_noedge,axis=1)
    profile3 = np.median(np.compress(spec1d_clip,spec2d_noedge,1),axis=1) # p in Dan's code

    xl = profile3.shape[0]
    x = np.arange(xl)

    guess,wstatus = param(profile3,x,xcen,verb)
    lsq = leastsq(gauss,guess,(profile3,x,verb),full_output=1,
                  xtol=1e-4,ftol=1e-4,gtol=1e-4,epsfcn=1e-4,factor=0.1)
    #              )


    params = lsq[0]
    status = lsq[4]

    if plot:
        fig = plt.figure()
        p = fig.add_subplot(111)
        #p.plot(spec2d_noedge)
        #p.plot(spec1d_noedge)
        #p.plot(spec1d_clip)
        #p.plot(profile1,drawstyle="steps")
        #p.plot(profile2,drawstyle="steps")
        p.plot(profile3,drawstyle="steps")

        gf = gauss(params,0,x,fitting=0)
        p.plot(gf)
        plt.show()

    if not (0 < params[0] < x[-1]): status = 0
    params[1] = np.power(10,params[1])          # amplitude
    params[2] = np.power(10,params[2])          # sigma
    fmt = "%12s %8.1f %4d" + 4 * "%9.1f" + " %r %r"
    out = (obj,D["flux"],xcen,) + tuple(params) + (status,wstatus,)
    print fmt % out
    D["info"] = params
    D["status"] = status
    if verb>1: raw_input("<< cr >> ")
    return params,status


    if not (0 < params[0] < x[-1]): status = 0
    params[1] = pow(10,params[1])
    params[2] = pow(10,params[2])
    fmt = "%12s %8.1f %4d" + 4 * "%9.1f" + " %r %r"
    out = (ob,D["flux"],h,) + tuple(params) + (status,wstatus,)
    print fmt % out
    D["info"] = params
    D["status"] = status
    if verb>1: raw_input("<< cr >> ")
    return params,status


def getdata(r,refobjs,ext=0,plot=0):

    objS = objD[r]["setup"]
    f = objS["spec2d"]
    nf = objS["spec2d_sig"]
    sf = objS["spec2d_sky"]
    ff = objS["spec2d_ext"]
    id = objS["object"]

    try:
        pfs = pyfits.open(f)
    except:
        print "Missing file %s!!!" % f
        D = {"id":id,"pos":None,"shift":True,"head":None,"data":None,
             "ref":None,"flux":None}
        return D

    spec2d = pfs[ext].data
    head = pfs[ext].header
    objold = head["objold"]
    #objlo  = head["objlo"]
    #objhi  = head["objhi"]
    objpos = head["objpos"]
    pfs.close()

    pfn = pyfits.open(nf)
    noise2d = pfn[ext].data
    pfn.close()

    pfs = pyfits.open(sf)
    sky2d = pfs[ext].data
    pfn.close()

    pff = pyfits.open(ff)
    flat2d = pff[ext].data
    pfn.close()

    # mini dictionary to store basic information
    D = {"id":id,"pos":objpos,"shift":True,"head":head}

    if D["id"] in refobjs: D["ref"] = True
    else: D["ref"] = False

    if 1 or D["ref"]:
        median_profile1 = np.median(spec2d,axis=1)
        median_profile2 = median_profile1 - np.percentile(median_profile1, 10)/2.0
        flux = np.sum(median_profile2)
        D["flux"] = flux
        #D["data"] = spec2d,noise2d
        N = spec2d.shape
        # between
        filt = np.where(np.greater(spec2d,-1e38)*np.less(spec2d,1e38))
        D["data"] = {}
        D["data"]["spec2d"] = spec2d[filt].reshape(N).astype("Float32")
        D["data"]["noise2d"] = noise2d[filt].reshape(N).astype("Float32")
        D["data"]["sky2d"] = sky2d[filt].reshape(N).astype("Float32")
        D["data"]["flat2d"] = flat2d[filt].reshape(N).astype("Float32")


        print "%15s %10.2f" % (id, flux)

        if plot:
            fig = plt.figure()
            p = fig.add_subplot(111)
            p.plot(median_profile1)
            p.plot(median_profile2)
            plt.show()

    return D

# rename to extraction
def Etrace(o,d,p,t,a,r=2,shift=True,f=False,edge=3):
    if not f: sys.stdout.write(" %r" % (o)); sys.stdout.flush()
    if f: s = d
    else: (s,n) = d
    x = arange(s.shape[1])
    y2,x2 = indices(s.shape,Float32)
    y2 = y2 - t
    if shift:
#       sm = VTKMedian(s,7,1,numret=1)

       ### dealing with severe cosmic rays
       sg = VTKGauss(s,sigmax=4,sigmay=1,numret=1)
       sd = s-sg
       ss = 1.49*VTKMedianSmooth(abs(sd),16,4,1)
       sr = divz(sd,ss)
       sb = 1.*equal(VTKConvolve(1.*less(sr,5),k=3,numret=1),1)
       sm = divz(VTKGauss(sb*s,sigmax=3,sigmay=1,numret=1),\
                       VTKGauss(sb,sigmax=3,sigmay=1,numret=1))

       for j in range(3):
         # typecode is weird, it thinks it is a byte, needs to be a float
         # this is essentially a mask
         u = 1.0*between(y2,p-max([8,a/2.0]),p+max([8,a/2.0]))
         uo = not_equal(sum(u),0)
         u = u * uo[NewAxis,::]
         #qdump("u.fits",u)
         #qdump("sp.fits",u*sm)
         #u = u * between(x2, s.shape[1]/4., 3*s.shape[1]/4.)
         sp = add.reduce(u*sm)
         np = sqrt(add.reduce(u*power(n,2)))
         sn = divz(sp,np)
         #print sn.shape
         sn = compress(sn,sn)
         #print sn.shape
         if len(sn) > len(s)/3: sn = getpct(sn,75)
         else: sn = 0.0
         #print
         #print sn,r
         if sn > r:
           um = u*not_equal(s,0)*not_equal(n,0)
           #qdump("um1.fits",um)
           ed = greater(y2,edge)*less(y2,y2.shape[0]-(1+edge))
           um = um*sm*greater(sm,-1*n)*ed
           #qdump("um2.fits",um)
           p2 = divz(sum((sm*um*y2).flat), sum((sm*um).flat))
           p2 = clip(p2,p-max([8,a]),p+max([8,a]))
           sys.stdout.write("(%.1f,%.1f,%.1f)" % (sn,p,p2-p)); sys.stdout.flush()
           p = p2
           #raw_input("<< cr >>")
         else:
           sys.stdout.write("(%.0f,%1.f,XXX)" % (sn,p)); sys.stdout.flush()
    spec = zeros(x.shape,Float64)
    noise = zeros(x.shape,Float64)
    y2 = y2 - p
    hi = greater(y2,-a/2.0)
    hi2 = greater(y2,-a/2.0-1) - hi
    hip = -(y2+a/2.0) * hi2
    lo = less(y2,a/2.0)
    lo2 = less(y2,a/2.0+1) - lo
    lop = (y2-a/2.0) * lo2
    mid = lo*hi
    use = (mid+hip+lop)
    intspec = lambda d: add.reduce(d*use)

    if f: return (intspec(s),p)
    else: return (intspec(s),getsqrt(intspec(power(n,2))),p)


##################################################
# usage:
#    extract_1dspec.py -smf 021953.SMF -dir 021953
##################################################

parser = argparse.ArgumentParser(description='Extract 1d spectra')

parser.add_argument('-smf', metavar='file', type=str, nargs='?',
                    default=None,help='SMF file')
parser.add_argument('-dir', metavar='path/to/dir', type=str, nargs='?',
                    default=None,help='directory of 2d spectra')
#parser.add_argument('-mark', metavar='file', type=bool, nargs='?',
#                    default=None,help='manually mark the object positions')
parser.add_argument('--trace-order', metavar='order of trace', type=int, nargs='?',
                    default=2,help='order of trace of atmospheric dispersion (global)')
parser.add_argument('--edge', metavar='pixels', type=float, nargs='?',
                    default=6.0,help='edge of the slit to exclude position fits [pixels]')
parser.add_argument('--trace-all', action="store_true",
                    help='trace all objects to for atmospheric dispersion (not reference objects)')
#parser.add_argument('--trace-ind', metavar='order of trace', type=int, nargs='?',
#                    default=2,help='trace each individual object (all) for atmospheric dispersion')
parser.add_argument('-o', metavar='order', type=int, nargs=2,
                    help="orders of fit to objects positions on the mask")
parser.add_argument('-a', metavar='aperture', type=int, nargs=1, default=10,
                    help="default aperture to extract 1D spectra [integer]")

args = parser.parse_args()

smf  = args.smf
fdir = args.dir
aper = args.a
edge = args.edge
orders = args.o
trace_all = args.trace_all
trace_order = args.trace_order

sw = 10

trace_ind = 0
trace_th = 2.0
flux_cut = 0.0
no_find = 0 # skip finding object positions
mark = 0

#sky = 1
#flat = 1

db = simpledb(smf,fdir)
objD = db.dict["objects"]
rows = db.row_search()
refobjs = db.refobjs
smf_ids = db.smf_ids
smf_x = db.smf_x
smf_y = db.smf_y

#print db
#print rows
#for d in objD: print d,objD[d]



D = [getdata(r,refobjs) for r in rows]
#print D

fluxes = np.array([d["flux"] for d in D if d["ref"]])
print fluxes

#if mark:
#else:
#[GetProfile(d,edge) for d in D if d["ref"] and d["flux"] > flux_cut]
[GetProfile(d,edge) for d in D if d["flux"] > flux_cut]

for d in D: d["ref"] = (d["flux"] > flux_cut)

sigmas = np.array([d["info"][2] for d in D if d["ref"]])
status = np.array([d["status"] for d in D if d["ref"]])
dpos = np.array([d["info"][0]-d["pos"] for d in D if d["ref"]])             # delta_pos?
print "median(dpos)=",np.median(dpos)

amps = np.array([d["info"][1] for d in D if d["ref"]])
pflux = np.array([d["flux"] for d in D if d["ref"]])                        # pflux?
refx = np.array([smf_x[smf_ids.index(d["id"])] for d in D if d["ref"]])
refy = np.array([smf_y[smf_ids.index(d["id"])] for d in D if d["ref"]])
refo = np.array([(d["id"] in refobjs) for d in D if d["ref"]])              # refobjs?

# PROBLEM!!!
# possible status error when no solutions are found
#good = np.greater(amps,0.0)*np.equal(status,1.0)
good = np.greater(amps,0.0)*np.not_equal(status,0)

print amps
print sigmas
# it is not doing what it was designed to do
#refo,refx,refy,pflux,amps,dpos,sigmas,status=np.compress(good,[refo,refx,refy,pflux,amps,dpos,sigmas,status],1)
print amps
print sigmas

ampcut = 0.0
weights = np.greater(amps,ampcut)
weights = weights*between(sigmas,0.1,2*np.median(sigmas))

print "made it"


# fit positions of the objects on the mask and correct the SMF file positions

# if the orders of the fit are already known
if orders: ox,oy = orders
else: ox,oy = 0,0
print ox,oy
print refx
print refy

fit_orders = True
#user_input = "go"
user_input = "done"

#    if not mark and (testsmf and not orders): user_input = "go"
#    else: user_input = "done"

#plt.ioff()
print plt.rcParams['interactive']

fig = plt.figure()
plt.show(block=False)

while fit_orders:

    markers = [['.','o'][r] for r in refo]
    colors = [['k','None'][r] for r in refo]
    colors_b = [['b','None'][r] for r in refo]
    colors_r = [['r','None'][r] for r in refo]

    p1 = fig.add_subplot(231)
    for i in xrange(len(refo)): p1.scatter(pflux[i],amps[i],c=colors[i],s=20, marker=markers[i],edgecolor='k')
    p1.set_xlabel("Flux")
    p1.set_ylabel("Amp")
    
    p2 = fig.add_subplot(232)
    p2.set_xlabel("X pos")
    p2.set_ylabel("$\sigma$")
    #p2.scatter(refx,sigmas,c="k",s=10)
    for i in xrange(len(refo)): p2.scatter(refx[i],sigmas[i],c=colors[i],s=20, marker=markers[i],edgecolor='k')
    p2.set_xlim(-100,100)
    p2.set_ylim(0.0,5.0)
    
    p3 = fig.add_subplot(233)
    p3.set_xlabel("X pos")
    p3.set_ylabel("Y pos")
    magnify = 10.0
    [p3.arrow(refx[j],refy[j],0,dpos[j]*magnify) for j in range(len(refx))] # need to fix arrors
    p3.set_xlim(-100,100)
    p3.set_ylim(-100,100)
    
    p4 = fig.add_subplot(234)
    p4.set_xlabel("X pos")
    p4.set_ylabel("$\Delta$Y pos")
    #p4.scatter(refx,dpos,c="k",s=10)
    for i in xrange(len(refo)): p4.scatter(refx[i],dpos[i],c=colors[i],s=20, marker=markers[i],edgecolor='k')
    p4.set_xlim(-100,100)
    p4.set_ylim(-20,20)
    
    lrefx = refx/300.0   # mm?
    lrefy = refy/300.0
    lxb = legendre(lrefx,ox)
    lyb = legendre(lrefy,oy)
    basis_fns = [xb*yb for xb in lxb for yb in lyb]
    
    print basis_fns
    
    uweights = 1*weights
    
    nfit = 1
    
    for fitpass in range(nfit):
       #sol = linear_least_squares(np.transpose(basis_fns)*uweights[::,np.newaxis],dpos*uweights)[0]
       sol = linalg.lstsq(np.transpose(basis_fns)*uweights[::,np.newaxis],dpos*uweights)[0]
       if not fitpass:
          # using Dan's planefit routine based on M-estimators 
          sol = planefit(np.compress(uweights,transpose(basis_fns),0),np.compress(uweights,dpos),sol,minfunc=mad)[0]
       if type(sol) is FloatType: sol = np.array([sol])
       fitted = np.dot(sol,basis_fns)
       resid = dpos - fitted
       if fitpass < nfit-1:
         uweights = uweights*np.less(np.abs(resid)/bwt(np.compress(uweights,resid))[1],5.0)
       #scatter = np.sqrt(np.sum(np.power(resid,2)*uweights)/np.sum(uweights))
       scatter = np.sqrt(divz(np.sum(np.power(resid,2)*uweights),np.sum(uweights)))
       print "SOL=",sol
       print "RMS=",scatter
       print "BWT=",bwt(compress(uweights,resid))
    # not sure why red is not ploting like imacsSpec1d.py 
    for i in xrange(len(refo)): p4.scatter(refx[i],resid[i],c=colors_r[i], s=20, marker=markers[i],edgecolor="r")

    refo_filt = np.compress(uweights,refo)
    refx_filt = np.compress(uweights,refx)
    resid_filt = np.compress(uweights,resid)
    for i in xrange(len(refo_filt)): p4.scatter(refx_filt[i],resid_filt[i],c=colors_b[i],s=20,marker=markers[i],edgecolor="b")
    #p4.scatter(refx,resid,c="r",s=10)
    #p4.scatter(np.compress(uweights,refx),np.compress(uweights,resid),c="b",s=10)
    
    scatter = np.sqrt(divz(np.sum(np.power(resid,2)*uweights),np.sum(uweights)))
    print "SOL=",sol
    print "RMS=",scatter
    print "BWT=",bwt(compress(uweights,resid))
    
    p4_label_r0 = "(ox,oy,cut)= (%d,%d,%f)" % (ox,oy,ampcut)
    p4_label_r1 = "SOL = %s" % (sol)
    p4_label_r2 = "RMS = %.6f" % (scatter)
    
    p4.text(0.05,0.95, p4_label_r0, transform=p4.transAxes)
    p4.text(0.05,0.90, p4_label_r1, transform=p4.transAxes)
    p4.text(0.05,0.85, p4_label_r2, transform=p4.transAxes)
    
    p5 = fig.add_subplot(235)
    p5.set_xlabel("Y pos")
    p5.set_ylabel("$\Delta$Y pos")
    #p5.scatter(refy,resid,c="k",s=10)
    for i in xrange(len(refo)): p5.scatter(refy[i],resid[i],c=colors[i],s=20, marker=markers[i],edgecolor='k')
    p5.set_xlim(-100,100)
    p5.set_ylim(-20,20)
    
    p6 = fig.add_subplot(236)
    p6.set_xlabel("X")
    p6.set_ylabel("Y")
    [p6.arrow(refx[j],refy[j],0,resid[j]*magnify) for j in range(len(refx))] # need to fix arrors
    p6.set_xlim(-100,100)
    p6.set_ylim(-100,100)

    plt.draw()
    
    if user_input != "done":
        prompt = "<(ox,oy,cut) = (%d,%d,%f)> " % (ox,oy,ampcut)
        user_input = raw_input( prompt)
        try:
           ox,oy,ampcut = eval(user_input)
           weights = np.greater(amps,ampcut)
           weights = weights*between(sigmas,0.1,2*np.median(sigmas))
        except: pass

    plt.clf()

    if user_input == "done":
      fig.savefig("%s_1dspec_leastsq.ps" % (smf))
      fit_orders = False
    
# Correcting postions based fit
lrefx = smf_x/300.0
lrefy = smf_y/300.0
lxb = legendre(lrefx,ox)
lyb = legendre(lrefy,oy)
basis_fns = [xb*yb for xb in lxb for yb in lyb]

dpos = dot(sol,basis_fns)
for d in D:
    if d["pos"] != None:
        j = smf_ids.index(d["id"])
        # original positions
        d["orig_pos"] = 1*d["pos"]
        if not mark and not no_find:
            # new positions
            d["pos"] = d["pos"] + dpos[j]

# Trace atmospheric dispersion
for d in D:
   if d["ref"] and (trace_all or d["id"] in refobjs) and (d["info"][1] > 0):
      d["info"][0] = d["pos"]
      d["trace"] = GetTrace(d,d["info"],edge,trace_order=trace_order,sw=sw)
   else:
      if trace_ind: d["trace"] = np.zeros(trace_order+1)
      else: d["trace"] = None

traces = np.asarray([d["trace"] for d in D if d["trace"] != None or trace_ind])
print traces

# Find the trace for each object
if trace_ind:
    trace = np.asarray([d["trace"] for d in D])
    if no_find:
        for tr in trace: tr[0] = 0.0
    print trace
# Find the global trace
else:
    trace = np.mean(traces,axis=0)
    trace[0] = 0.0
    print trace
    trace = [trace]

print
for j in range(len(traces)):
   print "Trace #%d = %r" % (j,traces[j])
if not trace_ind: print "Trace=",trace
print

#msigma = add.reduce(uweights*sigmas)/add.reduce(uweights)
msigma = divz(np.add.reduce(uweights*sigmas),np.add.reduce(uweights))
#print sigmas
print "msigma=",msigma
print

#sys.exit()

#if flat:
#    print "Reading %s" % (flat)
#    f = FITS(flat)
#    FD = [getextractdata(s,refobjs,flat) for s in 1+arange(f["N_SLITS"])]
#    f.close()
#    print "Done reading %s" % (flat)
#    for fd in FD:
#        if fd["ref"]: fd["ap"] = a
#        else: fd["ap"] = a
#        j = smf_ids.index(fd["id"])
#        if not no_find: fd["pos"] = fd["pos"] + dpos[j]


#if sky:
#    print "Reading %s" % (lamp)
#    f = FITS(lamp)
#    LD = [getextractdata(s,refobjs,lamp) for s in 1+arange(f["N_SLITS"])]
#    f.close()
#    print "Done reading %s" % (lamp)
#    for ld in LD:
#        if ld["ref"]: ld["ap"] = a
#        else: ld["ap"] = a
#        j = smf_ids.index(ld["id"])
#        if not no_find: ld["pos"] = ld["pos"] + dpos[j]


# Spectral Extractions
for d in D:
    if d["ref"]: d["ap"] = aper
    else: d["ap"] = aper




# pickle traces
naxis1 = D[0]["data"]["spec2d"].shape[1]
x0, x1 = 0.0, 1.0*naxis1
xc, xw = (x1+x0)/2.0, (x1-x0)/2.0
tfun = []
for tr in trace:
    b = legendre((np.arange(0,naxis1,1,Float32)-xc)/xw,len(tr)-1)
    tfun.append(dot(tr,b))

tfun = np.asarray(tfun)
#db.dict["trace"] = tfun
db.write()

#print tfun.shape
t = open(fdir + ".trace","wb")
cPickle.dump(tfun,t) 
t.close()



#
        #print D
        #print D[0].keys()
        #print len(D)

#        t0 = time.time()

        ###
        ### way of doing different size apertures for extracting data ###
        ###
    
#        if unique:
#            items = getunique(unique)
#            for item in items:
#                id,ap,pos = item
#                for d in D:
#                    if d["id"] == id:
#                        print "found one!"
#                        d["ap"] = ap
#                        d["pos"] = pos
#                        d["shift"] = False

 
print "Starting extractions"
spectra = [Etrace(d["id"],(d["data"]["spec2d"],d["data"]["noise2d"]),d["pos"],tfun[[0,k][trace_ind]],d["ap"],r=trace_th,shift=not no_find and d["shift"],f=False,edge=edge) for k,d in enumerate(D) if d["data"] != None]
print
# spectra, noise

# trace positions?
Sp = array([s[2] for s in spectra])
#print Sp
#print [1+d["orig_pos"] for d in D]
#print [1+d["pos"] for d in D]

spectra_arr = np.array([[s[0] for s in spectra],[s[1] for s in spectra]])

print spectra_arr.shape
#sys.exit()

# Sky Extractions
print "Starting sky extractions"
skys = [Etrace(d["id"],d["data"]["sky2d"],d["pos"],tfun[[0,k][trace_ind]],d["ap"],r=0.0,shift=False,f=True,edge=edge) for k,d in enumerate(D) if d["data"] != None]
skys_arr = np.array([s[0] for s in skys])

## Flattened-Flatfield Extractions
print "Starting flat extractions"
flats = [Etrace(d["id"],d["data"]["flat2d"],d["pos"],tfun[[0,k][trace_ind]],d["ap"],r=0.0,shift=False,f=True,edge=edge) for k,d in enumerate(D) if d["data"] != None]
flats_arr = np.array([f[0] for f in flats])

print flats_arr.shape
print spectra_arr.shape

fluxed_spectra_arr = np.array([[divz(s[0]/(f[0]/np.median(f[0]))) for s,f in zip(spectra,flats)],
                               [divz(s[1]/(f[0]/np.median(f[0]))) for s,f in zip(spectra,flats)]])

print fluxed_spectra_arr.shape

all_extractions = concatenate([spectra_arr,[skys_arr],[flats_arr],fluxed_spectra_arr],0)
print all_extractions.shape

#print D
#print len(D)


# individual 1dspec FITS files
k = 0
for d in D:
    if d["data"] != None:
        id = d["id"]
        head = d["head"]
        spectrum = all_extractions[::,k,::]
        
        nf = fdir + "/" + id + "_1dspec.fits"
        
        new_head = head.copy()
        #new_head["extrold"] = 1. + Sp[k] - d["orig_pos"]
        #new_head["extrpos"] = 1. + Sp[k] - d["pos"]
        #new_head["extrold"] = 1. + Sp[k] - d["orig_pos"]
        new_head["extrpos"] = (Sp[k],"1D extracted position of spectrum, in pixels")
        new_head["aper"] = (d["ap"],"aperture size, diameter in pixels")
        
        new_head["ARRAY1"] = ("SPECTRUM","units of counts")
        new_head["ARRAY2"] = ("NOISE","units of counts")
        new_head["ARRAY3"] = ("SKY","units of counts")
        new_head["ARRAY4"] = ("RAW FLATS","units of counts")
        new_head["ARRAY5"] = ("FLUXED SPECTRUM","units of counts")
        new_head["ARRAY6"] = ("FLUXED NOISE","units of counts")
        
        hdu = pyfits.PrimaryHDU(spectrum,header=new_head)
        #hdu = pyfits.PrimaryHDU(spectrum)
        hdu.writeto(nf,clobber=True)
        k += 1

#        
#        #########################
#        #########################
#        #########################
#        t1 = time.time()
#        print "Done with all extractions in %.3fs" % (t1-t0)
#        print
#        
#        naxes = array(c.shape)[::-1].tolist()
#
#        if testbig:
#            f = FITS(fitsfile,readonly=0)
#            apnums = [f["APNUM%d" % (1+j)].split() for j in range(len(D))]
#            apcens = array([f["cntrl%03d" % (1+j)] + (1 + Sp[j] - D[j]["pos"]) for j in range(len(D))])
##            apcens = array([f["csect%da" % (1+j)] + Sp[j] for j in range(len(D))])
#            aplo,aphi = apcens[NewAxis,::] + 0.5*a*array([-1,1])[::,NewAxis]
#            for j in range(len(D)): apnums[j][2:4] = "%.2f" % (aplo[j]), "%.2f" % (aphi[j])
#            fixpos = lambda j: f["cntrl%03d" % (1+j)] + (1 + Sp[j] - D[j]["orig_pos"])
##            fixpos = lambda j: f["csect%da" % (j+1)] + Sp[j]
#            [f.modkey("excen%03d" % (1+j), fixpos(j)) for j in range(len(D))]
#            f.close()
#
#        if testbig: ffitsfile = fitsfile.replace("_big.fits","_2p_1spec.fits")
#        else: ffitsfile = fitsfile.replace("_2spec.fits","_1spec.fits")

nf = fdir + "_1dspec.fits"
#hdu = pyfits.PrimaryHDU(S,header=head)
hdu = pyfits.PrimaryHDU(all_extractions)
hdu.writeto(nf,clobber=True)

#
#        ff=FITS(ffitsfile,readonly=0,create=1)
#        f = FITS(fitsfile)
#        if not testbig: f.goto(1)
#        ff.copyheader(f)
#        f.close()
#        ff.close()
#
#        ff=FITS(ffitsfile,readonly=0)
#        ff.resize(naxes)
#        #[ff.modkey(key,val) for (key,val) in [("dapert",a)]]
#        if testbig:
#          [ff.modkey("APNUM%d" % (1+j)," ".join(apnums[j])) for j in range(len(D))]
#          [ff.modkey("EXAP%03d" % (1+j),D[j]["ap"]) for j in range(len(D))]
#
#        if flat: ff.modkey("FLATFILE",flat,comment="")
#        if lamp: ff.modkey("LAMPFILE",lamp,comment="")
#        ff.flush()
#        ff.putdata(c)
#        ff.close()
#        print ffitsfile
#
#


























#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#    L_all = open(smffile,"r").readlines()
#    L = [l.split() for l in L_all]
#
#    instrument = [l[1] for l in L if l[0] == "INSTRUMENT"][0]
#
#    refobjs = [l[1] for l in L if l[0] == "HOLE"]
#    smf_ids = [l[1] for l in L if l[0] in ["HOLE","SLIT"]]
#    smf_x = array([float(l[4]) for l in L if l[0] in ["HOLE","SLIT"]])
#    smf_y = array([float(l[5]) for l in L if l[0] in ["HOLE","SLIT"]])
#
#    if not refobjs and not trace_all and not testsmf:
#        raise "No reference objects in list! You may need to use --trace-all"
#
#    t0 = time.time()
#    f = FITS(fitsfile)
#    if testbig:
#        D = [getbigdata(s,smf_ids,testsmf) for s in 1+arange(f["N_SLITS"])]
#    else:
#        D = [getdata(s,smf_ids,testsmf,instrument) for s in 1+arange(f["N_SLITS"])]

#    fluxes = array([d["flux"] for d in D if d["ref"]])




















###################################################################################################################
#
# >>> from basics import *
# >>> import inspect
# >>> inspect.getmodule(linear_least_squares)
# <module 'oldnumeric.linear_algebra' from '/data/software/CarPy/dist/lib/oldnumeric/linear_algebra.pyc'>
# >>> 
#
###################################################################################################################
