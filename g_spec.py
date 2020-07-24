#!/usr/bin/env python

# display 1d and 2d spectra

import glob,re,os,sys
from types import *
import argparse
import cPickle

from astropy.convolution import convolve, Box1DKernel, Box2DKernel
#import pyfits
import astropy.io.fits as pyfits
import numpy as np
import matplotlib
#print matplotlib.__version__
from scipy import integrate


# works!
matplotlib.use('Qt5Agg')
from PyQt5.QtCore import pyqtRemoveInputHook
pyqtRemoveInputHook()
# fixes the following:
#  QCoreApplication::exec: The event loop is already running
# during raw_input


# Frameworks builds?
#matplotlib.use('WXAgg')
#matplotlib.use('WX')

# doesn't work
#matplotlib.use('macosx')

# now fails
#matplotlib.use('TkAgg') # segmentation fault/ crashes gui/ broken now?

import matplotlib.pyplot as plt

#from matplotlib.widgets import TextBox

from matplotlib import patches
#plt.switch_backend('macosx')  # seems to work better with the mac
import matplotlib.cm as cm
from matplotlib.widgets import TextBox


#from image_utils import rebin2D_sum
#from math_utils import MAD
#from algorithms import rect_smooth, tri_smooth

from simpledb import simpledb

def rebin2D_sum(a, shape):
    #print shape
    #print a.shape[0]
    #print a.shape[0]//shape[0]
    # 2D
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).sum(-1).sum(1)

# from functons
def poly_n(p, data, x):
    p = np.array(p)
    if type(x) == FloatType: y = 0
    else: y = np.zeros((x.shape[0],))
    for n in np.arange(p.shape[0]):
        y += p[n]*x**n
    return y - data

def poly_x(x,p):
    return poly_n(p, 0, x)

def rect_smooth(y):
    j = 1
    go = 1
    S = np.zeros(y.shape)
    while go:
        S[j] = (y[j-1]+y[j]+y[j+1])/3
        j += 1
        if j == len(y) - 1: go = 0
    return S

def tri_smooth(y):
    j = 2
    go = 1
    S = np.zeros(y.shape)
    #print S.shape
    while go:
        S[j] = (y[j-2]+2*y[j-1]+3*y[j]+2*y[j+1]+y[j+2])/9
        j += 1
        if j == len(y) - 2: go = 0
    return S

def MAD(x):
    return np.nanmedian(np.abs(x-np.nanmedian(x)))


class SingleObject:

    def __init__(self,fig,fr,N,trace=None):
        self.fr = fr
        self.fig = fig

        self.fit = False
        self.fit_flux = []
        self.fit_wave = []
        self.fit_window = []
        self.fit_text = ""
        self.fit_output_fmt = "table"
        #self.fit_output_fmt = "single"

        self.slits = N
        
        #self.slits = header['NSLITS']

        self.redshift = 0
        self.z = None
        self.zq = None
        self.shift_is_held = False
        self.control_is_held = True
        self.sum_mode = "flux"
        self.s_is_held = False
        self.g_is_held = False
        self.window = []
        self.smooth = 1
        self.comment = ""
        self.verb = 0
        self.mask = 1

        self.trace = 0
        self.tfun = None

        self.spec_i = 0
        self.noise_i = 1

        self.display()

    def display(self,resetbounds = 1):

        #print plt.isinteractive()

        if resetbounds:
            self.setup_obj = objD[rows[self.fr]]["setup"]
            #self.smf_obj = objD[rows[self.fr]]["smf"]

            #self.tfun = objD["trace"] 
            self.tfun = db.dict["trace"]

            self.z = self.setup_obj["z"] 
            self.zq = self.setup_obj["zq"] 

            if self.setup_obj["spec1d_flux"]:
                f1d = self.setup_obj["spec1d_flux"]
                self.spec_i = 6
                self.noise_i = 7
            else:
                f1d = self.setup_obj["spec1d"]
                self.spec_i = 0
                self.noise_i = 1

            f2d = self.setup_obj["spec2d"]

            pf1d = pyfits.open(f1d,ignore_missing_end=1)
            
            self.spec1d = pf1d[0].data
            header = pf1d[0].header
            
            self.crpix1  = header['CRPIX1']     # starting pixel
            self.crval1  = header['CRVAL1']     # starting wavelength
            self.cd1_1   = header['CD1_1']      # dispersion
            self.dc_flag = header['DC-FLAG']    # Log-linear flag
            #print self.spec1d.shape
            
            ysize,xsize = self.spec1d.shape
            self.x = np.arange(xsize)+1
            self.w = (self.x-self.crpix1)*self.cd1_1+self.crval1  # compute wavelength
            
            #print (xsize-self.crpix1)*self.cd1_1+self.crval1
            if self.dc_flag: self.w = np.power(10,self.w)
            #print self.w
            
            #print self.spec1d.shape
            #print self.spec1d[self.spec_i,:].shape
            #print self.spec1d[self.spec_i,self.fr,:].shape
            #print self.spec1d[self.spec_i,self.fr,:]
            #print
            #print

            #self.extrold = header["extrold"] 
            self.objpos = header["objpos"] 
            self.extrpos = header["extrpos"]
            self.aper = header["aper"]
            if self.extrpos == "N/A": self.extrpos = self.extrold
            #print self.objpos,self.extrpos,self.aper
            
            pf2d = pyfits.open(f2d,ignore_missing_end=1)
            self.big = pf2d[0].data
            self.header = pf2d[0].header
            #print self.big.shape

            #self.ya    = self.header['CSECT%0dA' % (self.fr+1)]
            #self.yb    = self.header['CSECT%0dB' % (self.fr+1)]
            #self.obj   = self.header['OBJ%03d' % (self.fr+1)]
            #self.apnum = self.header['APNUM%d' % (self.fr+1)]

            self.obj = self.setup_obj["object"]

            #self.img   = self.big[self.ya:self.yb,:]
            self.img   = self.big
            #print self.spec1d.shape
            self.spec  = self.spec1d[self.spec_i,:]
            self.noise = self.spec1d[self.noise_i,:]

            rect_spec = rect_smooth(self.spec)
            tri_spec = tri_smooth(self.spec)
            self.smooth_spec = tri_spec
            self.smooth_img = self.img

            self.fit = False
            self.fit_flux = []
            self.fit_wave = []
            self.fit_window = []
            self.fit_text = ""

            self.bounds()


        #self.z = self.setup_obj["z"]

        #self.w1 = 10010

        self.fig.clf()
        #plt.clf()
        # [x0, y0, xwidth, ywidth]
        self.ax1 = plt.axes([0.1, 0.7, 0.8, 0.2])
        self.ax2 = plt.axes([0.1, 0.1, 0.8, 0.6], sharex=self.ax1)


        # working with making text box
        #axbox = plt.axes([0.1, -0.1, 0.8, 0.075])
        #initial_text=""
        #text_box = TextBox(axbox, 'Evaluate', initial=initial_text)
        #text_box.on_submit(submit)


        # link xaxis together for zoom

        #print self.big[self.ya:self.yb,:].shape

        self.x0 = int((self.w0 - self.crval1)/self.cd1_1+self.crpix1)
        self.x1 = int((self.w1 - self.crval1)/self.cd1_1+self.crpix1)


        #print self.x0,self.x1
        #print self.w0,self.w1
        #print self.img.shape
        #print self.spec.shape

        #self.ax1.imshow(self.img[:,self.x0:self.x1], origin='lower',
        #                interpolation="nearest", cmap=cm.gray_r,
        #                vmin=self.v0, vmax=self.v1, aspect="auto",
        #                extent=(self.w0,self.w1,0,self.img.shape[0]))
        self.ax1.imshow(self.smooth_img[:,self.x0:self.x1], origin='lower',
                        interpolation="nearest", cmap=cm.gray_r,
                        vmin=self.v0, vmax=self.v1, aspect="auto",
                        extent=(self.w0,self.w1,0,self.img.shape[0]))


        plt.setp(self.ax1.get_xticklabels(), visible=False)




        self.ax2.plot(self.w,self.spec,color="k",alpha=0.2,drawstyle="steps")
        #self.ax2.plot(self.w,tri_spec,color="k",drawstyle="steps")
       # self.ax2.plot(self.w,self.rect_spec,color="k",drawstyle="steps")
        self.ax2.plot(self.w,self.noise,color="r",drawstyle="steps")

        if self.smooth:
            self.ax2.plot(self.w,self.smooth_spec,c="k",drawstyle="steps-mid")

        #self.ax2.plot(self.x,self.spec1d[self.spec_i,self.fr,:],color="k",drawstyle="steps")
        #self.ax2.plot(self.x,self.spec1d[self.noise_i,self.fr,:],color="r",drawstyle="steps")

        self.ax1.set_xlim(self.w0,self.w1)
        self.ax2.set_xlim(self.w0,self.w1)


        self.ax2.set_ylim(self.y0,self.y1)

        self.ax2.set_xlabel("Observed Wavelength ($\AA$)",fontsize=12)
        if self.setup_obj["spec1d_flux"]:
            self.ax2.set_ylabel("Flux (erg/s/cm$^2$/$\AA$)",fontsize=12)
        else:
            self.ax2.set_ylabel("Flux (counts)",fontsize=12)


        # display info
        self.ax1.text(0.0,1.1,"row = %i" % rows[self.fr],
                      transform=self.ax1.transAxes)
        self.ax1.text(0.2,1.1,"id = %s" % self.obj,
                      transform=self.ax1.transAxes)
        if self.redshift and type(self.z) is not NoneType:
            self.display_redshift()
            self.ax1.text(0.70,1.1,"z = %.4f" % self.z,
                          transform=self.ax1.transAxes)
            self.ax1.text(0.90,1.1,"zq = %r" % self.zq,
                          transform=self.ax1.transAxes)

        if self.trace and self.tfun is not NoneType:
            self.display_trace()
            #self.ax1.text(0.85,1.1,"%.2f" % self.cntr_old,
            #    transform=self.ax1.transAxes)
            #self.ax1.text(0.85,1.1,"%.2f" % self.cntr_new,
            #    transform=self.ax1.transAxes)

        #print self.v0,self.v1

        if self.mask:
            self.mask_atmos()

        if self.fit:
            self.display_fit()


        # text box prototype (work in progress)
        #self.comment = objD[rows[self.fr]]["setup"]["comment"]
        #self.text_box = TextBox(ax2, 'Comment', initial=self.comment)
        #self.text_box.on_submit(submit)
        #print plt.isinteractive()

        #plt.ion()
        #plt.draw()

        self.fig.canvas.draw()
        #print plt.isinteractive()

    def display_trace(self):

        self.ax1.plot(self.w, self.extrpos + self.tfun[0],"--",color="w")
        self.ax1.plot(self.w, self.extrpos + self.tfun[0]-self.aper/2.0,"--",color="w")
        self.ax1.plot(self.w, self.extrpos + self.tfun[0]+self.aper/2.0,"--",color="w")
        #self.ax1.plot(self.w, self.extrold + self.tfun[0],"--",color="y")
        self.ax1.plot(self.w, self.objpos + self.tfun[0],"--",color="b")


    def display_fit(self):
        self.ax2.plot(self.fit_wave,self.fit_flux,c="g")

        #print self.fit_window
        x_ends = np.array(self.fit_window)[:,0]
        y_ends = np.array(self.fit_window)[:,1]

        self.ax2.scatter(x_ends,y_ends,c="g")
        self.ax2.text(0.02,0.8,self.fit_text, transform=self.ax2.transAxes,color="g")

    def display_redshift(self):
        dw = (self.w1-self.w0)/200.
        if len(wlist) > 0:
            for wl in wlist:
                wemit,n,line = wl
                wemit = float(wemit)
                n = int(n)
                c = ["r","b"][n-1]
                wobs = wemit*(1+self.z)
                if wobs > self.w0 and wobs < self.w1:
                    self.ax2.plot([wobs,wobs],[self.y0,self.y1],"--",color=c)
                    self.ax2.text(wobs-dw,0.85*self.y1,line,fontsize=fontsize,
                        rotation='vertical',
                        horizontalalignment='center',
                        verticalalignment='center',)
        else:
            print "Line list not loaded!"

#axbox = plt.axes([0.1, 0.05, 0.8, 0.075])
#text_box = TextBox(axbox, 'Evaluate', initial=initial_text)
#text_box.on_submit(submit)


    def mask_atmos(self):
        Aband = [7580,7700]    
        Bband = [6850,6950]
        Atmos = [Aband, Bband]

        for band in Atmos:
            #print band

            w0,w1 = band
        
            yh = self.y1-self.y0
            xw = w1-w0
            
            #xc = (w0+w1)/2.
            #yc = (y0+y1)/2.
        
            #print xc,yc
            #print xw,yh
        
            self.ax2.add_patch(patches.Rectangle(
                                           (w0, self.y0),   # (x,y)
                                            xw,          # width
                                            yh,          # height
                                            facecolor = "0.9",
                                            edgecolor = "0.9",
                                            zorder=5,
                                            alpha=0.85,
                                         )
                       )
         
        
    def bounds(self):
        sigma = 5

        vstd = 1.4826*MAD(self.big.flat)
        vmed = np.nanmedian(self.big.flat)
        ystd = 1.4826*MAD(self.spec1d[self.spec_i,:])
        ymed = np.nanmedian(self.spec1d[self.spec_i,:])

        #print ystd,ymed

        self.w0 = self.w[0]
        self.w1 = self.w[-1]

        self.x0 = self.x[0]
        self.x1 = self.x[-1]
 
        self.v0 = vmed-sigma*vstd
        self.v1 = vmed+sigma*vstd
        self.y0 = ymed-sigma*ystd
        self.y1 = ymed+sigma*ystd

    def on_key(self,event):
        print('you pressed', event.key, event.xdata, event.ydata)



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
        # test key events
        if self.verb:
            self.on_key(event)

        resetbounds = 0
        zf = 2.0
        
        #print self.v0,self.v1,self.w0,self.w1,self.y0,self.y1

        #if event.key == 'shift':
        #    self.shift_is_held = True
        #    print 'shift!'

        #if event.key == 'control':
        #    self.control_is_held = True
        #    print 'control!'

        #if event.key == "n":
        #   if self.shift_is_held:
        #       if self.fr > 0:
        #           self.fr -= 1
        #           resetbounds = 1
        #   else:
        #       if self.fr < self.slits-1:
        #           self.fr += 1
        #           resetbounds = 1

        if event.key == 'g':
            val = int(raw_input("Go to frame? "))
            if val >= 0 and val < self.slits:
                self.fr = val
                resetbounds = 1
            else:
                print "Error! Out of range!"

        if event.key == "n":
           if self.fr < self.slits-1:
               self.fr += 1
               resetbounds = 1

        if event.key == "N":
           if self.fr > 0:
               self.fr -= 1
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

        # Boxcar smoothing
        # http://joseph-long.com/writing/AstroPy-boxcar/
        if event.key == 'B':
            val = int(raw_input("Enter the boxcar width? "))
            self.smooth_spec = convolve(self.spec, Box1DKernel(val))
            self.smooth_img = convolve(self.img,  Box2DKernel(val))

        if event.key == 'R':
            binx,biny = map(int, raw_input("Enter the rebinning (x,y)? ").split(',') )
            print binx,biny

            yl,xl = self.img.shape

            if xl % binx: nx = xl/binx
            else: nx = xl/binx

            if yl % biny: ny = yl/biny
            else: ny = yl/biny
            #print

            #print nx,ny
            

            self.smooth_img = rebin2D_sum(self.img[:ny*biny,:nx*binx],(ny,nx))


        if event.key == 'T':
            if self.trace: self.trace = 0
            else: self.trace = 1

        #if event.key == "x" or event.key == "z":
        #    if not self.control_is_held:
        #      if self.shift_is_held:
        #          xran = xran/zf; self.w0 = xc - xran/2.; self.w1 = xc + xran/2.
        #          if self.w0 < self.w[0]:  self.w0 = self.w[0]
        #          if self.w1 > self.w[-1]: self.w1 = self.w[-1]
        #      else:
        #          xran = xran*zf; self.w0 = xc - xran/2.; self.w1 = xc + xran/2.
        #          if self.w0 < self.w[0]:  self.w0 = self.w[0]
        #          if self.w1 > self.w[-1]: self.w1 = self.w[-1]

        #if event.key == "y" or event.key == "z":
        #    if not self.control_is_held:
        #      if self.shift_is_held:
        #          yran = yran/zf
        #          self.y0 = yc - yran/2.; self.y1 = yc + yran/2.
        #      else:
        #          yran = yran*zf
        #          self.y0 = yc - yran/2.; self.y1 = yc + yran/2.

        if event.key == "h":
            self.w0 = self.w0 - xran/10.
            if self.w0 < self.w[0]: self.w0 = self.w[0]
            self.w1 = self.w0 + xran
        if event.key == "l":
            self.w1 = self.w1 + xran/10.
            if self.w1 > self.w[-1]: self.w1 = self.w[-1]
            self.w0 = self.w1 - xran
        if event.key == "j":
            self.y0 = self.y0 - yran/10.
            self.y1 = self.y1 - yran/10.
        if event.key == "k":
            self.y0 = self.y0 + yran/10.
            self.y1 = self.y1 + yran/10.

        if event.key == "r":
            resetbounds = 1

        if event.key == "w":
            try:
                wav_obs = float(raw_input("Enter rest-frame wavelength of line? "))
                self.z = event.xdata/wav_obs - 1
                objD[rows[self.fr]]["setup"]["z"] = self.z
                #self.resetbounds = 1
            except:
                print "Error!"

        if event.key == "z":
            if self.redshift: self.redshift = 0
            else: self.redshift = 1

        if event.key == "H":
            zq_th = 1
            bins = 20

            zs = np.array([objD[i]["setup"]["z"] for i in objD])
            zqs = np.array([objD[i]["setup"]["zq"] for i in objD])

            filt = (zs != np.array(None))*(zqs != np.array(None))*((zqs > zq_th))

            new_z = zs[filt]

            fig2 = plt.figure()
            fig2.canvas.set_window_title(fdir)
            p2 = fig2.add_subplot(111)
            p2.hist(new_z,bins)
            p2.text(0.0,1.025,"N = %i" % (len(new_z)), transform=p2.transAxes)
            p2.text(0.9,1.025,"zq > %i " % (zq_th), transform=p2.transAxes)
            p2.set_xlabel("Redshift z")

            #plt.draw()
            plt.show()


        if event.key == "M":
           print "Mark new position"
           print objD[rows[self.fr]]
           print xc, yc
           #os.system('extract_1dspec.py -dir %s -smf %s -obj %s -pos %s' % (fdir,smf,obj,pos))
           print

        if event.key == "A":
           print "Adjust aperture"
           val = int(raw_input("Enter new aperture (FWHM) in pixels? "))
           print objD[rows[self.fr]]
           print xc, yc
           #os.system('extract_1dspec.py -dir %s -smf %s -obj %s -ap %s' % (fdir,smf,obj,ap))
           print

        if event.key == "E":
           choice = "Extract additional object at %.1f? (y/n)" % (yc)

           if choice == "y":
               aper = int(raw_input("Enter new aperture (FWHM) in pixels? "))
               print objD[rows[self.fr]]
               print xc, yc
               #os.system('extract_1dspec.py -dir %s -smf %s -newobj %s' % (fdir,smf,obj))
               print


   
        if event.key == "Z":
            try:
                self.z = float(raw_input("Enter redshift? "))
                objD[rows[self.fr]]["setup"]["z"] = self.z
                #self.resetbounds = 1
            except:
                print "Error!"

        if event.key == "Q":
            print "Quiting..."
            db.write()
            sys.exit()

        #if event.key in ["[","]"] and wlist and self.redshift:
        if event.key in "[];',.{}" and wlist and self.redshift:
            self.z = self.setup_obj["z"]
            if type(self.z) is NoneType: self.z = 0

            if event.key == "[": self.z -=0.1
            if event.key == "]": self.z +=0.1
            if event.key == ";": self.z -=0.01
            if event.key == "'": self.z +=0.01
            if event.key == ",": self.z -=0.001
            if event.key == ".": self.z +=0.001
            if event.key == "{": self.z -=0.0001
            if event.key == "}": self.z +=0.0001

            objD[rows[self.fr]]["setup"]["z"] = self.z
            
        #if event.key in ["[","]"] and wlist and self.redshift:
        #    self.z = self.setup_obj["z"]
        #    if type(self.z) is NoneType: self.z = 0

        #    if event.key == "[": self.z -=0.01
        #    if event.key == "]": self.z +=0.01

        #    objD[rows[self.fr]]["setup"]["z"] = self.z

        #if event.key in ["{","}"] and wlist and self.redshift:
        #    self.z = self.setup_obj["z"]
        #    if type(self.z) is NoneType: self.z = 0

        #    if event.key == "{": self.z -=0.001
        #    if event.key == "}": self.z +=0.001

        #    objD[rows[self.fr]]["setup"]["z"] = self.z

        # Redshift Quality
        #if event.key == "0":
        #    self.zq = 0
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq
        #if event.key == "1":
        #    self.zq = 1
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq
        #if event.key == "2":
        #    self.zq = 2
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq
        #if event.key == "3":
        #    self.zq = 3
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq
        #if event.key == "4":
        #    self.zq = 4
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq
        #if event.key == "5":
        #    self.zq = 5
        #    objD[rows[self.fr]]["setup"]["zq"] = self.zq

        if event.key == "q":
            try:
                self.zq = raw_input("Enter z-quality? ")
                objD[rows[self.fr]]["setup"]["zq"] = int(self.zq)
            except:
                print "Error!"

        # Predefined Redshifts
        if event.key == ")":
            self.z = 0.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "!":
            self.z = 1.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "@":
            self.z = 2.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "#":
            self.z = 3.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "$":
            self.z = 4.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "%":
            self.z = 5.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "^":
            self.z = 6.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "&":
            self.z = 7.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "*":
            self.z = 8.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "(":
            self.z = 9.00
            objD[rows[self.fr]]["setup"]["z"] = self.z


        if event.key == "C":
            try:
                self.comment = raw_input("Enter comment? ")
                objD[rows[self.fr]]["setup"]["comment"] = self.comment 
            except:
                print "Error!"

        #if event.key == "S":
        #    try:
        #        self.spectrum = raw_input("Classify spectrum? (e.g. C, E, C+A, C+E, C+A+E) ")
        #        objD[rows[self.fr]]["setup"]["spectrum"] = self.spectrum
        #    except:
        #        print "Error!"

        if event.key == 'S':
            # sum flux between two points
            # subtract off fit continuum
            if self.s_is_held:
                self.window.append([xc,yc])
                self.window.sort()
     
                self.fit_window = self.window

                print "(x,y) =",self.fit_window
    
                w0,y0 = self.window[0]
                w1,y1 = self.window[1]
    
                # integer data positions
                x0 = np.argmin(np.abs(self.w - w0))
                x1 = np.argmin(np.abs(self.w - w1))
    
                w0 = self.w[x0]
                w1 = self.w[x1]
     
                #print w0,w1
                #print x0,x1
                        
    
                # determine coefficients from two points
                m = (y1-y0)/(w1-w0)
                b = y0 - m*w0
                param = [b,m]

                wfit  = self.w[x0:x1+1]
                dw = wfit[1]-wfit[0]

                #print w0,w1
                #print y0,y1
                #print param
                #print dw
                #print
    
                # sum flux
                if self.sum_mode == "flux":
                    xfit_flux  = self.spec[x0:x1+1]
                    print 
                    #print integrate.trapz(xfit_flux,wfit)  
                    all_flux = np.sum(xfit_flux,0)*dw
                    cont_flux = integrate.quad(poly_x,w0,w1,args=(param))[0]
                    final_flux = all_flux - cont_flux
                    #print xfit_flux

                    x_ends = np.array(self.fit_window)[:,0]
                    y_ends = np.array(self.fit_window)[:,1]


                    xfit_err = self.noise[x0:x1+1]**2
                    #xfit_err = self.spec[x0:x1+1]**2
                    print np.sqrt(integrate.trapz(xfit_err,wfit))
                    all_err = np.sqrt(np.sum(xfit_err,0)*dw)
                    print all_err

 
                    if self.fit_output_fmt == "table":
                        #print "F(all)    F(cont)   F(line)"
                        #print "%8.2e  %8.2e  %8.2e" % (all,cont,final_flux)
                        print "%3s  %6s  %7s  %7s  %7s  %7s  %8s  %8s  %8s" % ("Row","ID","Wav1","Wav2","F1","F2","F(all)","F(cont)","F(line)")
                        print "%3i  %6s  %7.1f  %7.1f  %7.1f  %7.1f  %8.2e  %8.2e  %8.2e+-%8.2e" % (rows[self.fr],self.obj,x_ends[0],x_ends[1],y_ends[0],y_ends[1],all_flux,cont_flux,final_flux,all_err)
                    elif self.fit_output_fmt == "single":
                        print "Flux (all)  =", all
                        print "Flux (cont) =", cont
                        print "Flux (line) =", final_flux
                        print 
                    
                    # plot fits
                    pfit = poly_n(param,0,wfit)  # poly fit
                    #self.ax2.plot(wfit,pfit,c="g")
                    #self.fig.canvas.draw()

                    self.fit = True
                    self.fit_wave = wfit
                    self.fit_flux = pfit
                    self.fit_text = "F(all) = %.2e\nF(cont) = %.2e\nF(line) = %.2e$\pm$%.2e" % (all_flux,cont_flux,final_flux,all_err)
                    #print wfit
                    #print pfit




                # rms mode
                if self.sum_mode == "rms":
                    xfit  = self.spec[x0:x1+1]**2
                    print np.sqrt(integrate.trapz(xfit,wfit))
                    all   = np.sqrt(np.sum(xfit,0)*dw)
                    print all
    
                self.s_is_held = False
                self.window = []
            else:
                self.window = [[xc,yc]]
                #print self.window
                print "Press shift+s again"
                self.s_is_held = True

        if event.key == 'p':

            # cols = ["ra","dec","object","type","z","zq","flag","features", "comment"]

            cols = ["ra","dec","object","z","zq","spectrum","comment"]
            fmt  = ["%s","%s","%-18s","%7.4f","%6i","%-6s","%s"]

            for i in objD:

                # test output
                # print objD[i]["setup"]

                print "%3i" % i,

                for k,j in enumerate(cols):
                    try:
                        val = objD[i]["setup"][j]
                        print fmt[k] % (val),
                    except:
                        print " "*6,
                    #except KeyError:
                    #    print " "*6,
                print

        if event.key == 'X':

            reg_f = "redshifts.reg"
            print "writing %s" % (reg_f)
            f = open(reg_f,"w")
            f.write("fk5\n")
            for i in objD:
                ra = objD[i]["setup"]["ra"]
                dec = objD[i]["setup"]["dec"]
                z = objD[i]["setup"]["z"]
                zq = objD[i]["setup"]["zq"]
                obj = objD[i]["setup"]["object"]
                comment = objD[i]["setup"]["comment"]

                if z: z_str = "%.3f" % z
                else: z_str = "" 
                if comment: comment_str = "%s" % comment
                else: comment_str = "" 

                colors = ["green","yellow","red"]
                if zq > len(colors)-1 or zq < 0: zq = 0
                color = colors[zq]


                f.write('circle(%s,%s,2") # text={%s} color=%s\n' % (ra,dec,z_str,color))
                #f.write('circle(%s,%s,3") # text={%s}\n' % (ra,dec,comment_str))
                f.write('circle(%s,%s,3") # text={%s}\n' % (ra,dec,obj))
            f.close

            #os.system('xpaset -p ds9 regions load all %s &' % reg_f)


        self.display(resetbounds=resetbounds)

    #def on_button_press(self,event):
    #    #print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % ('double' if event.dblclick else 'single', event.button, event.x, event.y, event.xdata, event.ydata))
    #    resetbounds = 0
    #    
    #    if event.dblclick:
    #        try:
    #            wav_obs = float(raw_input("Enter rest-frame wavelength of line? "))
    #            self.z = event.xdata/wav_obs - 1
    #            objD[rows[self.fr]]["setup"]["z"] = self.z
    #            #self.resetbounds = 1
    #        except:
    #            print "Error!"
    #    
    #    self.display(resetbounds=resetbounds)


    def on_key_release(self, event):
        if event.key == 'shift':
            self.shift_is_held = False
        if event.key == 'control':
            self.control_is_held = False

    def connect(self):
        # http://matplotlib.org/users/event_handling.html
        'connect to all the events we need'
        #self.cid_key_press      = self.fig.canvas.mpl_connect(
        #     'button_press_event', self.on_button_press)
        self.cid_key_press      = self.fig.canvas.mpl_connect(
             'key_press_event', self.on_key_press)
        #self.cid_button_release = self.fig.canvas.mpl_connect(
        #     'button_release_event', self.on_button_release)
        self.cid_key_release    = self.fig.canvas.mpl_connect(
             'key_release_event', self.on_key_release)

    #def on_button_release(self,event):
    #    resetbounds = 0
    #    self.w0,self.w1 = self.ax1.get_xlim()                        
    #    self.y0,self.y1 = self.ax1.get_ylim()
    #    self.display(resetbounds=resetbounds)



# e.g
#  g_spec.py -smf 021953.SMF -dir 021953
#  g_spec.py -smf 021953.SMF -dir 021953 --trace 021953.trace

#try:
#    ph1 = os.getenv("PYTHONHOME1")
#    if not ph1: raise "PYTHONHOME1 environmental variable not set!"
#    wavefile = "/".join([ph1,"/datafiles/linelists/galaxylines.dat"])
#except:
#    wavefile = None

for p in sys.path:
    wavefile = p + "/data/galaxylines.dat"
    #wavefile = p + "/data/galaxylines_short.dat"
    #wavefile = p + "/../datafiles/linelists/galaxylines.dat"
    #wavefile = p + "/../datafiles/linelists/galaxylines_short.dat"
    if os.path.exists(wavefile): break
    else: wavefile = None


# raw
#spec_i = 0
#noise_i = 1
# "fluxed" divided by raw flat
#spec_i = 4
#noise_i = 5
# flux calibrated
#spec_i = 6
#noise_i = 7


# comparison
#spec_i = 0
#noise_i = 4 # fluxed spectra



# old format
#spec_i = 5
#noise_i = 6

mode = "manual"

fontsize=12

parser = argparse.ArgumentParser(description='Display 1d and 2d spectra.')

#parser.add_argument('-1d', metavar='file', type=str, nargs='?',
#                    help='1d spectra')
#parser.add_argument('-2d', metavar='file', type=str, nargs='?',
#                    help='2d spectra')
parser.add_argument('prefix', metavar='file', type=str, nargs='?',
                    default=None,help='prefix for SMF/save/trace files')
parser.add_argument('-smf', metavar='file', type=str, nargs='?',
                    default=None,help='SMF file')
parser.add_argument('-fr', metavar='frame', type=int, nargs='?',
                    default=0, help='frame number')
parser.add_argument('-l', metavar='file', type=str, nargs='?',
                    default=wavefile, help='line list')
parser.add_argument('-dir', metavar='file', type=str, nargs='?',
                    default=None,help='directory of 2d spectra')
parser.add_argument('--trace', metavar='file', type=str, nargs='?',
                    default=None,help='trace file')

#args = parser.parse_args(namespace=)
args = parser.parse_args()
#print args
#print dir(args)
#print args.__dict__.keys()

prefix = args.prefix

if prefix != None:
    smf = prefix + ".SMF"
    fdir = prefix
    trace = prefix + ".trace"
else:
    smf = args.smf
    fdir = args.dir
    if smf == None and fdir == None:
       print "Need at least one argument!"
       sys.exit()

    trace = args.trace


#print prefix


fr = args.fr




db = simpledb(smf,fdir)
print db
objD = db.dict["objects"]
rows = db.row_search()
print rows
for d in objD: print d,objD[d]
print objD[0]['setup'].keys()

N = len(rows)


### linelists 
if wavefile and os.path.exists(wavefile):
    wlines = open(wavefile,"r").readlines()
    wlines = filter(lambda l: l[0] != "#", wlines)
    wlist = map(str.split, wlines)
else:
    wlist = []


# import traces
tfun = None
ind = 0
if trace:
    if trace[-6:] == ".trace":
        tfun = cPickle.load(open(trace))
        # Compatability mode (older version)
        #if type(tfun) is ListType: ind = 0
        # Temporary mode, future mode to be added to fs
        #if type(tfun) is np.ArrayType:
        #    if tfun.shape[0] > 1: ind = 1
        #    elif tfun.shape[0] == 1: ind = 0
        #print type(tfun)
        #print tfun.shape

    elif trace[-5:] == ".fits":
        pf = pyfits.open(trace,ignore_missing_end=1)
        tfun = [pf[0].data]
        #print trace[-5:]
        #print tfun

    else:
        print "trace format not supported"

#objD["trace"] = tfun
#print objD["trace"]
db.dict["trace"] = tfun





# reset matplotlib standard keymap
for i in plt.rcParams:
    if "keymap" in i:
        plt.rcParams[i] = ''

#sys.exit()

#plt.ion()

fig = plt.figure()
fig.canvas.set_window_title('g_spec.py')
#fig = plt.figure(figsize=(10,8))
SO = SingleObject(fig,fr,N)
SO.connect()
#fig.canvas.show()
plt.show()
#plt.draw()




# convert to FITS
objD = db.dict["objects"]
rows = db.row_search()
#print rows

row = []
id = []
oclass = []
redshift = []
quality = []
comment = []
extpos = []
extaper = []
extflag = []
alignbox = []


print
print
print
#print
#print

for i,d in enumerate(objD): 

    row.append(i+1)
    id.append(objD[d]["setup"]["object"])

    if objD[d]["setup"]["smf_type"] == "SLIT": 
        oclass.append("galaxy")
    elif objD[d]["setup"]["smf_type"] == "HOLE": 
        oclass.append("star")

    if objD[d]["setup"]["z"] == None:
        redshift.append(0)
    else:
        redshift.append(objD[d]["setup"]["z"])

    if objD[d]["setup"]["zq"] == None:
        quality.append(0)
    else:
        quality.append(objD[d]["setup"]["zq"])


    comment.append(objD[d]["setup"]["comment"])
 
    #print d,objD[d]

    f1d = objD[d]["setup"]["spec1d"]
    if f1d == None:
        extpos.append(-99)
        extaper.append(-99)
    else:
        pf1d = pyfits.open(f1d,ignore_missing_end=1)
        head = pf1d[0].header
    
        extpos.append(head["EXTRPOS"])
        extaper.append(head["APER"])

    extflag.append(False)
    alignbox.append(False)

#print objD[0]['setup'].keys()


#print row
#print id 
#print oclass
#print redshift
#print quality
#print comment
#print extpos
#print extaper
#print extflag
#print alignbox

fits_output = fdir.replace("/","") + "_objects.fits"

# FITS table format
#row id class redshift quality comment extpos extaper extflag alignbox
col1 = pyfits.Column(name='row', format='K', array=np.array(row))
col2 = pyfits.Column(name='id', format='20A', array=id)
col3 = pyfits.Column(name='class', format='6A', array=oclass)
col4 = pyfits.Column(name='redshift', format='D', array=np.array(redshift))
col5 = pyfits.Column(name='quality', format='K', array=np.array(quality))
col6 = pyfits.Column(name='comment', format='100A', array=comment)
col7 = pyfits.Column(name='extpos', format='D', array=np.array(extpos))
col8 = pyfits.Column(name='extaper', format='D', array=np.array(extaper))
col9 = pyfits.Column(name='extflag', format='L', array=extflag)
col10 = pyfits.Column(name='alignbox', format='L', array=alignbox)

cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])

hdu = pyfits.BinTableHDU.from_columns(cols)
#hdu = pyfits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
hdu.writeto(fits_output,clobber=True)









db.write()
