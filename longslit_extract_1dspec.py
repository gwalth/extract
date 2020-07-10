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
import numpy.linalg as linalg

from scipy import signal
from scipy.optimize import leastsq
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from simpledb import simpledb
from math_utils import MAD

# https://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy
#def consecutive(data, stepsize=1):
#    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

#consecutive(a)
class GUI_mark:
    def __init__(self,fig,spec2d,profile):
        self.fig = fig
        self.spec2d = spec2d
        self.profile = profile

        self.crpix1,self.crval1,self.cd1_1 = wcs

        xsize = self.profile.shape[0]
        zsize,ysize = self.spec2d.shape
        self.zsize = zsize
        self.x = np.arange(xsize)+1

        print zsize,ysize

        self.shift_is_held = False
        self.control_is_held = False

        self.display()

    def display(self,resetbounds = 1):

        if resetbounds:
            self.bounds()


        #self.fig.clf()
        plt.clf()
        # [x0, y0, xwidth, ywidth]
        self.ax1 = plt.axes([0.1, 0.1, 0.2, 0.8])
        self.ax2 = plt.axes([0.3, 0.1, 0.6, 0.8])

        #self.x0 = 0.
        #self.x1 = 


        print self.x0,self.x1
        print self.w0,self.w1

        self.ax1.imshow(self.spec2d[:,self.x0:self.x1],
                        interpolation="nearest", cmap=cm.gray_r,
                        vmin=self.v0, vmax=self.v1, aspect="auto",
                        extent=(self.w0,self.w1,self.zsize,0))
        #self.ax1.imshow(self.spec2d[:,self.x0:self.x1],
        #self.ax1.imshow(self.spec2d,
        #                interpolation="nearest", cmap=cm.gray_r)
       
        # example code: 
        #   ax1.imshow(spec2d[:,x0:x1], interpolation="nearest", cmap=cm.gray_r, vmin=v0, vmax=v1, aspect="auto", extent=(w0,w1,zsize,0))

        plt.setp(self.ax1.get_xticklabels(), visible=False)

        self.ax1.set_xlim(self.w0,self.w1)
        self.ax2.set_xlim(self.w0,self.w1)



        self.ax2.plot(self.wav,self.spec1d,color="k",alpha=0.2,drawstyle="steps")


        self.ax2.set_ylim(self.y0,self.y1)

        self.ax2.set_xlabel("Observed Wavelength ($\AA$)",fontsize=12)
        self.ax2.set_ylabel("Flux",fontsize=12)


        # display info
        #self.ax1.text(0.0,1.1,"id = %s" % self.obj,
        #              transform=self.ax1.transAxes)
    def bounds(self):
        sigma = 5

        # filter out zeros
        fwav = np.compress(np.greater(self.wav,0),self.wav)

        self.x0 = np.argmin(np.abs(self.wav - self.w0))
        self.x1 = np.argmin(np.abs(self.wav - self.w1))

        vstd = 1.4826*MAD(self.spec2d.flat)
        vmed = np.median(self.spec2d.flat)

        ystd = 1.4826*MAD(self.profile[self.x0:self.x1])
        ymed = np.median(self.profile[self.x0:self.x1])
        #self.x0 = self.x[0]
        #self.x1 = self.x[-1]
 
        print self.x0
        print self.x1

        self.v0 = vmed-sigma*vstd
        self.v1 = vmed+sigma*vstd

        self.y0 = ymed-sigma*ystd
        self.y1 = ymed+sigma*ystd

        print self.v0
        print self.v1
        print self.y0
        print self.y1

    def on_key_press(self,event):

        xc, yc = event.xdata, event.ydata

        #plt.gcf()
        #print plt.gcf()
        #self.fig.gcf()

        #print self.obj
        #print self.fr
        #print self.ax1.get_xlim() 
        #print self.ax2.get_xlim() 
        self.w0,self.w1 = self.ax2.get_xlim() 
        self.y0,self.y1 = self.ax2.get_ylim() 

        #print event.key
        resetbounds = 0
        zf = 2.0
        
        #print self.v0,self.v1,self.w0,self.w1,self.y0,self.y1

        if event.key == 'shift':
            self.shift_is_held = True

        if event.key == 'control':
            self.control_is_held = True

        if event.key == "n":
           if self.shift_is_held:
               if self.fr > 0:
                   self.fr -= 1
                   resetbounds = 1
           else:
               if self.fr < self.slits-1:
                   self.fr += 1
                   resetbounds = 1

        vran = self.v1 - self.v0
        yran = self.y1 - self.y0
        xran = self.w1 - self.w0

        if event.key == "s": # sharper
           self.v0 = self.v0 + vran/4.
           self.v1 = self.v1 - vran/4.
        if event.key == "d": # duller
           self.v0 = self.v0 - vran/4.
           self.v1 = self.v1 + vran/4.
        if event.key == "b": # brighter
           self.v0 = self.v0 + vran/4.
           self.v1 = self.v1 + vran/4.
        if event.key == "f": # fainter
           self.v0 = self.v0 - vran/4.
           self.v1 = self.v1 - vran/4.
        #if event.key == "i": # invert 
        #   self.v0,self.v1 = self.v1,self.v0

        if event.key == "x" or event.key == "z":
            if not self.control_is_held:
              if self.shift_is_held:
                  xran = xran/zf; self.w0 = xc - xran/2.; self.w1 = xc + xran/2.
                  if self.w0 < self.w[0]:  self.w0 = self.w[0]
                  if self.w1 > self.w[-1]: self.w1 = self.w[-1]
              else:
                  xran = xran*zf; self.w0 = xc - xran/2.; self.w1 = xc + xran/2.
                  if self.w0 < self.w[0]:  self.w0 = self.w[0]
                  if self.w1 > self.w[-1]: self.w1 = self.w[-1]

        if event.key == "y" or event.key == "z":
            if not self.control_is_held:
              if self.shift_is_held:
                  yran = yran/zf
                  self.y0 = yc - yran/2.; self.y1 = yc + yran/2.
              else:
                  yran = yran*zf
                  self.y0 = yc - yran/2.; self.y1 = yc + yran/2.

        if event.key == "h":
            self.x0 = self.x0 - xran/10.
            if self.x0 < self.x[0]: self.x0 = self.x[0]
            self.x1 = self.x0 + xran
        if event.key == "l":
            self.x1 = self.x1 + xran/10.
            if self.w1 > self.w[-1]: self.x1 = self.x[-1]
            self.x0 = self.x1 - xran
        if event.key == "j":
            self.y0 = self.y0 - yran/10.
            self.y1 = self.y1 - yran/10.
        if event.key == "k":
            self.y0 = self.y0 + yran/10.
            self.y1 = self.y1 + yran/10.

        if event.key == "r":
            resetbounds = 1

        if event.key == "z":
            if self.control_is_held:
                if self.redshift: self.redshift = 0
                else: self.redshift = 1


        self.display(resetbounds=resetbounds)

    def on_key_release(self, event):
        if event.key == 'shift':
            self.shift_is_held = False
        if event.key == 'control':
            self.control_is_held = False

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
    #print x0,x1
    #print wt
    #print wt2
    #print kz

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
        #sys.exit()

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
    print j, xcen, amp
    if amp > 0: amp = np.log10(amp)
    else: amp = 2.0
    sig = np.log10(1.5)
    print sig
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

def GetTrace(D,params,edge,trace_order=2,sw=10,plot=1):
    #print D["data"]
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


    #naxis1 = D["data"]["spec2d"].shape[1]
    #x0, x1 = 0.0, 1.0*naxis1
    #xc, xw = (x1+x0)/2.0, (x1-x0)/2.0
    #b = legendre((np.arange(0,naxis1,1,Float32)-xc)/xw,len(traces)-1)



    if plot:
        #vstd = 1.4826*MAD(spec2d.flat)
        #vmed = np.median(spec2d.flat)

        #print np.percentile(spec2d,1)
        #print np.percentile(spec2d,99)
        #print np.min(spec2d)
        #print np.max(spec2d)
        #print np.median(spec2d.flat)


        per_lo = np.percentile(spec2d,0.5)
        per_hi = np.percentile(spec2d,99.5)

        #print np.sum(spec2d > per_hi)
        #print np.sum(spec2d < per_lo)

        #print "I am here!"
        fig = plt.figure()
        p = fig.add_subplot(111)

        p.imshow(spec2d, origin='lower',
                 interpolation="nearest", cmap=cm.gray_r,
        #         aspect="auto")
        #         aspect="auto",vmin=-2*vstd,vmax=5*vstd)
                 aspect="auto", vmin=per_lo, vmax=per_hi)

        xr = x[0,:]
        br = legendre((xr-xc)/xw,trace_order)
        tfun = dot(t,br)
        #print tfun
        #print tfun.shape
        #p.plot(xr, tfun+params[0],"--",color="r")

        p.plot(xr, tfun+params[0],"--",color="b",label="Extraction aperture")
        p.plot(xr, tfun+params[0]-aper/2.0,"--",color="b")
        p.plot(xr, tfun+params[0]+aper/2.0,"--",color="b")

        p.plot(xr, tfun+params[0]-large_aper/2.0,"--",color="c",label="Large extraction aperture")
        p.plot(xr, tfun+params[0]+large_aper/2.0,"--",color="c")
        p.set_xlabel("Y-axis (pixels)")
        p.set_ylabel("X-axis (pixels)")
        leg = p.legend(loc=2,numpoints=1)

        plt.show()

    return t
#    Test the solution by rebinning spectra:
#    p = reshape(dot(t,b),x.shape)
#    r = VTKImageTransform(spec2d,0.0*p,p,reverse=0,numret=1)

def GetProfile(D,edge,plot=0,verb=0,mark=0):
    #(spec2d,noise2d) = D["data"]
    spec2d = D["data"]["spec2d"]
    noise2d = D["data"]["noise2d"]
    N = spec2d.shape
    y, x = np.indices(N,dtype=np.float32) 
    print N
    EdgeMask = np.logical_or(np.less(y,edge),np.greater(y,y.shape[0]-(1+edge))) 
    spec2d_noedge = spec2d*np.logical_not(EdgeMask)                       # w in Dan's code

    xcen = D["pos"]
    obj = D["id"]
    spec1d_noedge = np.sum(spec2d_noedge,axis=0)                          # pz in Dan's code
    spec1d_clip = np.not_equal(spec1d_noedge,0)                           # gz in Dan's code

    #profile1 = np.median(spec2d,axis=1)
    #profile2 = np.median(spec2d_noedge,axis=1)
    profile3 = np.median(np.compress(spec1d_clip,spec2d_noedge,1),axis=1) # p in Dan's code
    profile4 = np.sum(np.compress(spec1d_clip,spec2d_noedge,1),axis=1) # p in Dan's code

    xl = profile3.shape[0]
    x = np.arange(xl)

    #guess,wstatus = param(profile3,x,xcen,verb=verb)
    guess,wstatus = param(profile4,x,xcen,verb=verb)
    print
    print guess
    #lsq = leastsq(gauss,guess,(profile3,x,verb),full_output=1,
    lsq = leastsq(gauss,guess,(profile4,x,verb),full_output=1,
                  xtol=1e-4,ftol=1e-4,gtol=1e-4,epsfcn=1e-4,factor=0.1)
    #              )


    params = lsq[0]
    status = lsq[4]
    print params
    print
    xcen = params[0]
    sigma = params[1]
    FWHM = 2.35*sigma
    print
    print "     pixels arcsec"
    print "FWHM  %.2f  %.2f" % (FWHM,FWHM*pix_scale)

    if plot:
        fig = plt.figure()
        p1 = fig.add_subplot(111)
        #p1.plot(spec2d_noedge)
        #p1.plot(spec1d_noedge)
        #p1.plot(spec1d_clip)
        #p1.plot(profile1,drawstyle="steps")
        #p1.plot(profile2,drawstyle="steps")
        #p1.plot(x,profile3,drawstyle='steps-mid')
        p1.plot(x,profile4,drawstyle='steps-mid',label="Profile of spectrum (sum)",c="k")

        dx = 0.1
        gx = np.arange(0,xl,dx)
        gf = gauss(params,0,gx,fitting=0)



        g1 = list(gx).index(int(xcen-aper/2.0+0.5))
        g2 = list(gx).index(int(xcen+aper/2.0+0.5))

        l1 = list(gx).index(int(xcen-large_aper/2.0+0.5))
        l2 = list(gx).index(int(xcen+large_aper/2.0+0.5))

        #print gx[g1:g2]
        #print gf[g1:g2]

        # total integrated 
        total_int = integrate.trapz(gf,gx)  
        #total_sum = np.sum(profile3,0)*1.0
        total_sum = np.sum(profile4,0)*1.0
        # large aperture integrated
        large_int = integrate.trapz(gf[l1:l2],gx[l1:l2])  
        #large_sum = np.sum(profile3[xcen-large_aper/2.0:xcen+large_aper/2.0],0)*1.0
        large_sum = np.sum(profile4[xcen-large_aper/2.0:xcen+large_aper/2.0],0)*1.0
        # normal aperture integrated
        aper_int = integrate.trapz(gf[g1:g2],gx[g1:g2])  
        #aper_sum = np.sum(profile3[xcen-aper/2.0:xcen+aper/2.0],0)*1.0
        aper_sum = np.sum(profile4[xcen-aper/2.0:xcen+aper/2.0],0)*1.0

        print
        print "       pixels arcsec"
        print "Aper_1    %3i  %4.1f" % (aper,pix_scale*aper)
        print "Aper_2    %3i  %4.1f" % (large_aper,pix_scale*large_aper)
        print
        print "       Aper_1    Aper_2    Aper_3"
        print "       [%i pix]  [%i pix]  [total]" % (aper,large_aper)
        print "F_sum  %.2e  %.2e  %.2e" % (aper_sum,large_sum,total_sum)
        print "F_int  %.2e  %.2e  %.2e" % (aper_int,large_int,total_int)
        print 
        print "                     sum   int"
        print "F(Aper_1)/F(Aper_2)  %.2f  %.2f" % (aper_sum/large_sum,aper_int/large_int) 
        print " [%i pix] [%i pix]" % (aper,large_aper)
        print "F(Aper_1)/F(Aper_3)  %.2f  %.2f" % (aper_sum/total_sum,aper_int/total_int)
        print " [%i pix] [total]" % (aper)
        print "F(Aper_2)/F(Aper_3)  %.2f  %.2f" % (large_sum/total_sum,large_int/total_int)
        print " [%i pix] [total]" % (large_aper)
        print
        print



        
        #print integrate.quad(poly_x,w0,w1,args=(param))[0]

        p1.plot(gx,gf,label="Gaussian fit")

        ymin,ymax = p1.get_ylim()

        p1.plot([xcen,xcen],[ymin,ymax],"--",c="b",label="Extraction aperture")
        p1.plot([xcen-aper/2.0,xcen-aper/2.0],[ymin,ymax],"--",c="b")
        p1.plot([xcen+aper/2.0,xcen+aper/2.0],[ymin,ymax],"--",c="b")
        p1.plot([xcen-large_aper/2.0,xcen-large_aper/2.0],[ymin,ymax],"--",c="c",label="Large extraction aperture")
        p1.plot([xcen+large_aper/2.0,xcen+large_aper/2.0],[ymin,ymax],"--",c="c")

        p1.set_ylim(ymin,ymax)
        p1.set_xlabel("X-axis (pixels)")
        p1.set_ylabel("Flux (counts)")
        leg = p1.legend(loc=2,numpoints=1)
        plt.show()

    if mark > 0:
        weights = divz(1,noise2d).flatten()**2
        F = spec2d.flatten()
        cf = divz(sum(F*weights),sum(weights))
        cf = 0.0
        rf = sqrt(divz(sum(weights*(F-cf)**2),sum(weights)))
        a1,a2 = cf-3*rf, +3*rf
        So = sort(F)
        No = len(So)
        a1,a2 = So[int(0.05*No+1)], So[int(0.95*No+1)]
        c1,c2 = 0,spec2d.shape[1]-1
        xfm = array([-0.5,1.0,0.0,-0.5,0.0,1.0])
        kcur = None
        xcur = spec2d.shape[1]/2
        ycur = params[0]
        #if testmark:
        #    for m in mf:
        #        if D["id"] == m[0].strip():
        #            params[0] = float(m[1])
        #            break
        #else:
        if 1:
            while kcur not in ["a","M","q"]:
                PG.cpgsvp(0.05,0.20,0.15,0.90)
                PG.newbox([0,spec2d.shape[0]-1],[min(p)-(max(p)-min(p))/10,max(p)+(max(p)-min(p))/3],ch=2.3,manvp=1)
                PG.xlabel("Y",ch=2.3)
                PG.ylabel("Profile",ch=2.3)
                PG.drawln(x,p,ls=1,hist=1)
                xa = arange(min(x),max(x)+1,0.01,Float)
                PG.drawln(x,p,ls=1,ci=2,hist=1)
                PG.drawln(xa,gauss(params,p,xa,0,0),ls=1,ci=4,hist=0)
                PG.drawln([params[0],params[0]],[min(p)-2*(max(p)-min(p)),max(p)+2*(max(p)-min(p))],ls=1,ci=7)
                PG.subtitle(D["id"],ch=2.3)

                PG.cpgsvp(0.25,0.95,0.15,0.90)
                PG.newbox([c1,c2],[0,spec2d.shape[0]-1],ch=2.3,manvp=1,era=0)
                PG.xlabel("X",ch=2.3)
                PG.ylabel("Y",ch=2.3)
                PG.subtitle(D["id"],ch=2.3)
                PG.cpggray(spec2d,spec2d.shape[1],spec2d.shape[0],\
                        1,spec2d.shape[1],1,spec2d.shape[0],a1,a2,xfm)
                PG.drawln([0,spec2d.shape[1]],[0.5+params[0],0.5+params[0]],ls=1,ci=7)
                PG.drawln([0,spec2d.shape[1]],[0.5+params[0]-10**params[2],0.5+params[0]-10**params[2]],ls=2,ci=7)
                PG.drawln([0,spec2d.shape[1]],[0.5+params[0]+10**params[2],0.5+params[0]+10**params[2]],ls=2,ci=7)
                PG.drawln([0,spec2d.shape[1]],[0.5+h,0.5+h],ls=1,ci=2)
                PG.drawln([0,spec2d.shape[1]],[0.5+h+mark/2.,0.5+h+mark/2.],ls=2,ci=2)
                PG.drawln([0,spec2d.shape[1]],[0.5+h-mark/2.,0.5+h-mark/2.],ls=2,ci=2)

                kcur,(xcur,ycur),(xdum,ydum) = PG.getcurs(xcur,ycur,xfm)
                if kcur == "I": raise "UserInterrupt"
                if kcur in ["a","q"]: break
                if kcur.upper() == "M":
                    params = array([ydum,99,3,0.0])
                    status = 1
                    if kcur == "M": break
              
                aran = a2-a1
                xran = c2-c1
                if kcur == "j": params[0] = max([0.0, params[0] - 0.5])
                if kcur == "k": params[0] = min([params[0] + 0.5, spec2d.shape[0]-1])
                if kcur == "J": params[0] = max([0.0, params[0] - 5.0])
                if kcur == "K": params[0] = min([params[0] + 5.0, spec2d.shape[0]-1])
                if kcur == "h":
                    c1 = c1 - xran/10.
                    if c1 < 0: c1 = 0
                    c2 = c1 + xran
                if kcur == "l":
                    c2 = c2 + xran/10.
                    if c2 > spec2d.shape[1]-1: c2 = spec2d.shape[1]-1
                    c1 = c2 - xran
              
                zf = 2.0
                if kcur == "s": a1=a1+aran/4.; a2=a2-aran/4.
                if kcur == "d": a1=a1-aran/4.; a2=a2+aran/4.
                if kcur == "f": a1=a1-aran/4.; a2=a2-aran/4.
                if kcur == "b": a1=a1+aran/4.; a2=a2+aran/4.
                if kcur == "r": a1,a2 = So[int(0.05*No+1)], So[int(0.95*No+1)]; c1,c2 = 0,spec2d.shape[1]-1
                if kcur == "x" or kcur == "X":
                    if kcur == "x": xran = xran*zf; c1 = xcur - xran/2.; c2 = xcur + xran/2.
                    if kcur == "X": xran = xran/zf; c1 = xcur - xran/2.; c2 = xcur + xran/2.
                    if c1 < 0: c1 = 0
                    if c2 > spec2d.shape[1]-1: c2 = spec2d.shape[1]-1
                PG.cpgeras()
            
            ### cheat we subtract a -1 to offset the -1 we would get in stage-extract-1dspec   
            params[0] -= 1.0
            mf.write("%-8s = %f\n" % (D["id"],params[0]))



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

def getdata_longslit(ext=0,plot=0):

    f = sys.argv[1]
    nf = sys.argv[2]
    sf = sys.argv[3]
    ff = sys.argv[4]

    id = "test"

    try:
        pfs = pyfits.open(f)
    except:
        print "Missing file %s!!!" % f
        D = {"id":id,"pos":None,"shift":True,"head":None,"data":None,
             "ref":None,"flux":None}
        return D

    spec2d = pfs[ext].data
    head = pfs[ext].header
    pfs.close()

    pfn = pyfits.open(nf)
    #noise2d = pfn[ext].data
    #print noise2d
    #noise2d = np.sqrt(pfn[ext].data)  # not exactly sure what the rsum frame is, i think thise should be close to the noise
    noise2d = getsqrt(pfn[ext].data)
    #print noise2d
    #sys.exit()
    pfn.close()

    pfs = pyfits.open(sf)
    sky2d = pfs[ext].data
    pfn.close()

    pff = pyfits.open(ff)
    flat2d = pff[ext].data
    fltexp = pff[ext].header['exptime']
    pfn.close()

    median_profile1 = np.median(spec2d,axis=1)
    median_profile2 = median_profile1 - np.percentile(median_profile1, 10)/2.0
    flux = np.sum(median_profile2)

    ###############################################################
    # Need rough position of object
    ###############################################################
    ###############################################################
    ###############################################################
    #objpos = 1000 # guess
    objpos = np.argmax(median_profile2)
    ###############################################################
    ###############################################################
    ###############################################################

    # mini dictionary to store basic information
    D = {"id":id,"pos":objpos,"shift":True,"head":head}
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

    D["fltexp"] = fltexp


    print "%15s %10.2f" % (id, flux)

    if plot:
        fig = plt.figure()
        p1 = fig.add_subplot(111)
        p1.plot(median_profile1,label="Profile of spectrum (median)",c="k")
        #p1.plot(median_profile2)
        p1.set_xlabel("X-axis (pixels)")
        p1.set_ylabel("Flux (counts)")
        leg = p1.legend(loc=2,numpoints=1)
        plt.show()

    return D

# rename to extraction
def Etrace(o,d,p,t,a,r=2,shift=True,f=False,edge=3,plot=0):
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

    if plot:
        fig = plt.figure()
        p1 = fig.add_subplot(111)
        p1.plot(intspec(s),drawstyle='steps-mid',label="Extracted spectrum",c="k")
        if not f: 
            p1.plot(getsqrt(intspec(power(n,2))),drawstyle='steps-mid',label="Extracted spectrum",c="r")
        #p1.set_xlabel("X-axis (pixels)")
        p1.set_ylabel("Flux (counts)")
        leg = p1.legend(loc=2,numpoints=1)
        plt.show()

    if f: return (intspec(s),p)
    else: return (intspec(s),getsqrt(intspec(power(n,2))),p)


##########################################################################################################
# usage:
#    #                          science              sigma              sky                  raw-flat
#    longslit_extract_1dspec.py ccd2051c_or1sum.fits ccd2051c_rsum.fits ccd2051c_mr1sum.fits ccd2055c_mrsum.fits
#    longslit_extract_1dspec.py ccd2052c_or1sum.fits ccd2052c_rsum.fits ccd2052c_mr1sum.fits ccd2053c_mrsum.fits
#    longslit_extract_1dspec.py ccd2091c_or1sum.fits ccd2091c_rsum.fits ccd2091c_mr1sum.fits ccd2099c_mrsum.fits

# longslit_extract_1dspec.py ccd2077c_or1sum.fits ccd2077c_rsum.fits ccd2077c_mr1sum.fits ccd2082c_mrsum.fits
##########################################################################################################

#parser = argparse.ArgumentParser(description='Extract 1d spectra')
#
#parser.add_argument('-smf', metavar='file', type=str, nargs='?',
#                    default=None,help='SMF file')
#parser.add_argument('-dir', metavar='path/to/dir', type=str, nargs='?',
#                    default=None,help='directory of 2d spectra')
##parser.add_argument('-mark', metavar='file', type=bool, nargs='?',
##                    default=None,help='manually mark the object positions')
#parser.add_argument('--trace-order', metavar='order of trace', type=int, nargs='?',
#                    default=2,help='order of trace of atmospheric dispersion (global)')
#parser.add_argument('--edge', metavar='pixels', type=float, nargs='?',
#                    default=6.0,help='edge of the slit to exclude position fits [pixels]')
#parser.add_argument('--trace-all', action="store_true",
#                    help='trace all objects to for atmospheric dispersion (not reference objects)')
##parser.add_argument('--trace-ind', metavar='order of trace', type=int, nargs='?',
##                    default=2,help='trace each individual object (all) for atmospheric dispersion')
#parser.add_argument('-o', metavar='order', type=int, nargs=2,
#                    help="orders of fit to objects positions on the mask")
#parser.add_argument('-a', metavar='aperture', type=int, nargs=1, default=10,
#                    help="default aperture to extract 1D spectra [integer]")

#args = parser.parse_args()

#smf  = args.smf
#fdir = args.dir
#aper = args.a
#edge = args.edge
#orders = args.o
#trace_all = args.trace_all
#trace_order = args.trace_order

sw = 10

trace_ind = 0
trace_th = 2.0
flux_cut = 0.0
no_find = 0 # skip finding object positions
mark = 0

edge = 6.0
trace_order = 2
#aper = 10
#aper = 12
aper = 14
#large_aper = 50
large_aper = 100

pix_scale = 0.189 # arcsec/pixel LDSS3

#sky = 1
flat = 1


# could make this work for multiple objects in a slit
D = getdata_longslit(plot=1)
#print D


#if mark:
#else:
#[GetProfile(d,edge) for d in D if d["ref"] and d["flux"] > flux_cut]
GetProfile(D,edge,plot=1)



sigmas = D["info"][2]
status = D["status"]
dpos = D["info"][0]
amps = D["info"][1]

# PROBLEM!!!
# possible status error when no solutions are found
#good = np.greater(amps,0.0)*np.equal(status,1.0)
good = np.greater(amps,0.0)*np.not_equal(status,0)


ampcut = 0.0
weights = np.greater(amps,ampcut)
weights = weights*between(sigmas,0.1,2*np.median(sigmas))

#print "made it"

# Trace atmospheric dispersion
D["trace"] = GetTrace(D,D["info"],edge,trace_order=trace_order,sw=sw)

traces = D["trace"]
#print traces


# Spectral Extractions
D["ap"] = aper




# pickle traces
naxis1 = D["data"]["spec2d"].shape[1]
x0, x1 = 0.0, 1.0*naxis1
xc, xw = (x1+x0)/2.0, (x1-x0)/2.0
tfun = []
b = legendre((np.arange(0,naxis1,1,Float32)-xc)/xw,len(traces)-1)
#print b
#print b.shape
tfun.append(dot(traces,b))

tfun = np.asarray(tfun)

id = D["id"]


#print tfun.shape
tf = open(id + ".trace","wb")
cPickle.dump(tfun,tf) 
tf.close()

 
print "Starting extractions"
spectra = [Etrace(D["id"],(D["data"]["spec2d"],D["data"]["noise2d"]),D["pos"],tfun[[0,0][trace_ind]],D["ap"],r=trace_th,shift=not no_find and D["shift"],f=False,edge=edge,plot=1)]
print

Sp = array([s[2] for s in spectra])

spectra_arr = np.array([[s[0] for s in spectra],[s[1] for s in spectra]])
print spectra_arr.shape


# Sky Extractions
print "Starting sky extractions"
#skys = [Etrace(d["id"],d["data"]["sky2d"],d["pos"],tfun[[0,k][trace_ind]],d["ap"],r=0.0,shift=False,f=True,edge=edge) for k,d in enumerate(D) if d["data"] != None]
skys = [Etrace(D["id"],D["data"]["sky2d"],D["pos"],tfun[[0,0][trace_ind]],D["ap"],r=0.0,shift=False,f=True,edge=edge)]
skys_arr = np.array([s[0] for s in skys])

## Flattened-Flatfield Extractions
if flat:
    print "Starting flat extractions"
    #flats = [Etrace(d["id"],d["data"]["flat2d"],d["pos"],tfun[[0,k][trace_ind]],d["ap"],r=0.0,shift=False,f=True,edge=edge) for k,d in enumerate(D) if d["data"] != None]
    flats = [Etrace(D["id"],D["data"]["flat2d"],D["pos"],tfun[[0,0][trace_ind]],D["ap"],r=0.0,shift=False,f=True,edge=edge)]
    flats_arr = np.array([f[0] for f in flats])

    print flats_arr.shape
    print spectra_arr.shape

    fluxed_spectra_arr = np.array([[divz(s[0]/(f[0]/np.median(f[0]))) for s,f in zip(spectra,flats)],
                               [divz(s[1]/(f[0]/np.median(f[0]))) for s,f in zip(spectra,flats)]])

    print fluxed_spectra_arr.shape

    all_extractions = concatenate([spectra_arr,[skys_arr],[flats_arr],fluxed_spectra_arr],0)

else:
    all_extractions = concatenate([spectra_arr,[skys_arr]],0)
print all_extractions.shape

#print D
#print len(D)


# individual 1dspec FITS files
id = D["id"]
head = D["head"]
spectrum = all_extractions[::,0,::]

nf = id + "_1dspec.fits"

new_head = head.copy()
#new_head["extrold"] = 1. + Sp[k] - d["orig_pos"]
#new_head["extrpos"] = 1. + Sp[k] - d["pos"]
#new_head["extrold"] = 1. + Sp[k] - d["orig_pos"]
new_head["extrpos"] = (Sp[0],"1D extracted position of spectrum, in pixels")
new_head["aper"] = (D["ap"],"aperture size, diameter in pixels")

new_head["ARRAY1"] = ("SPECTRUM","units of counts")
new_head["ARRAY2"] = ("NOISE","units of counts")
new_head["ARRAY3"] = ("SKY","units of counts")
if flat:
    new_head["ARRAY4"] = ("RAW FLATS","units of counts")
    new_head["ARRAY5"] = ("FLUXED SPECTRUM","units of counts")
    new_head["ARRAY6"] = ("FLUXED NOISE","units of counts")

    new_head["FLTEXP"] = (D["fltexp"],"flat exposure time")

hdu = pyfits.PrimaryHDU(spectrum,header=new_head)
#hdu = pyfits.PrimaryHDU(spectrum)
hdu.writeto(nf,clobber=True)


