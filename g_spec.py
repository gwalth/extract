#!/usr/bin/env python

# display 1d and 2d spectra

import glob,re,os,sys
from types import *
import argparse
import cPickle

from astropy.convolution import convolve, Box1DKernel
import pyfits
import numpy as np
import matplotlib
#print matplotlib.__version__
#matplotlib.use('TkAgg') # segmentation fault
#matplotlib.use('WXAgg')
#matplotlib.use('WX')
import matplotlib.pyplot as plt
#plt.switch_backend('macosx')  # seems to work better with the mac
import matplotlib.cm as cm


#from math_utils import MAD
#from algorithms import rect_smooth, tri_smooth

from simpledb import simpledb

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
    return np.median(np.abs(x-np.median(x)))


class SingleObject:

    def __init__(self,fig,fr,N,trace=None):
        self.fr = fr
        self.fig = fig


        self.slits = N
        
        #self.slits = header['NSLITS']

        self.redshift = 0
        self.z = None
        self.zq = None
        self.shift_is_held = False
        self.control_is_held = True
        self.smooth = 1
        self.comment = ""

        self.trace = 0
        self.tfun = None

        self.display()

    def display(self,resetbounds = 1):

        if resetbounds:
            self.setup_obj = objD[rows[self.fr]]["setup"]
            #self.smf_obj = objD[rows[self.fr]]["smf"]

            #self.tfun = objD["trace"] 
            self.tfun = db.dict["trace"]

            self.z = self.setup_obj["z"] 

            f1d = self.setup_obj["spec1d"]
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
            
            print self.spec1d.shape
            #print self.spec1d[spec_i,:].shape
            #print self.spec1d[spec_i,self.fr,:].shape
            #print self.spec1d[spec_i,self.fr,:]
            #print
            #print

            #self.extrold = header["extrold"] 
            self.objpos = header["objpos"] 
            self.extrpos = header["extrpos"]
            self.aper = header["aper"]
            if self.extrpos == "N/A": self.extrpos = self.extrold
            
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
            self.spec  = self.spec1d[spec_i,:]
            self.noise = self.spec1d[noise_i,:]

            rect_spec = rect_smooth(self.spec)
            tri_spec = tri_smooth(self.spec)
            self.smooth_spec = tri_spec

            self.bounds()


        #self.z = self.setup_obj["z"]

        #self.w1 = 10010

        #self.fig.clf()
        plt.clf()
        # [x0, y0, xwidth, ywidth]
        self.ax1 = plt.axes([0.1, 0.7, 0.8, 0.2])
        self.ax2 = plt.axes([0.1, 0.1, 0.8, 0.6], sharex=self.ax1)
        # link xaxis together for zoom

        #print self.big[self.ya:self.yb,:].shape

        self.x0 = int((self.w0 - self.crval1)/self.cd1_1+self.crpix1)
        self.x1 = int((self.w1 - self.crval1)/self.cd1_1+self.crpix1)


        #print self.x0,self.x1
        #print self.w0,self.w1
        #print self.img.shape
        #print self.spec.shape

        #self.ax1.imshow(self.big[self.ya:self.yb,self.x0:self.x1],
        self.ax1.imshow(self.img[:,self.x0:self.x1], origin='lower',
                        interpolation="nearest", cmap=cm.gray_r,
                        vmin=self.v0, vmax=self.v1, aspect="auto",
                        extent=(self.w0,self.w1,0,self.img.shape[0]))

        #self.ax1.imshow(self.big[self.ya:self.yb,self.x0:self.x1],
        #                interpolation="nearest", cmap=cm.gray_r,
        #                vmin=self.v0, vmax=self.v1, aspect="auto",
        #                extent=(self.x0,self.x1,self.ya,self.yb))

        plt.setp(self.ax1.get_xticklabels(), visible=False)




        self.ax2.plot(self.w,self.spec,color="k",alpha=0.2,drawstyle="steps")
        #self.ax2.plot(self.w,tri_spec,color="k",drawstyle="steps")
       # self.ax2.plot(self.w,self.rect_spec,color="k",drawstyle="steps")
        self.ax2.plot(self.w,self.noise,color="r",drawstyle="steps")

        if self.smooth:
            self.ax2.plot(self.w,self.smooth_spec,c="k",drawstyle="steps-mid")

        #self.ax2.plot(self.x,self.spec1d[spec_i,self.fr,:],color="k",drawstyle="steps")
        #self.ax2.plot(self.x,self.spec1d[noise_i,self.fr,:],color="r",drawstyle="steps")

        self.ax1.set_xlim(self.w0,self.w1)
        self.ax2.set_xlim(self.w0,self.w1)


        self.ax2.set_ylim(self.y0,self.y1)

        self.ax2.set_xlabel("Observed Wavelength ($\AA$)",fontsize=12)
        self.ax2.set_ylabel("Flux",fontsize=12)


        # display info
        self.ax1.text(0.0,1.1,"row = %i" % rows[self.fr],
                      transform=self.ax1.transAxes)
        self.ax1.text(0.2,1.1,"id = %s" % self.obj,
                      transform=self.ax1.transAxes)
        if self.redshift and type(self.z) is not NoneType:
            self.display_redshift()
            self.ax1.text(0.85,1.1,"z = %.4f" % self.z,
                          transform=self.ax1.transAxes)

        if self.trace and self.tfun is not NoneType:
            self.display_trace()
            #self.ax1.text(0.85,1.1,"%.2f" % self.cntr_old,
            #    transform=self.ax1.transAxes)
            #self.ax1.text(0.85,1.1,"%.2f" % self.cntr_new,
            #    transform=self.ax1.transAxes)

        #print self.v0,self.v1

        plt.draw()

    def display_trace(self):

        self.ax1.plot(self.w, self.extrpos + self.tfun[0],"--",color="w")
        self.ax1.plot(self.w, self.extrpos + self.tfun[0]-self.aper/2.0,"--",color="w")
        self.ax1.plot(self.w, self.extrpos + self.tfun[0]+self.aper/2.0,"--",color="w")
        #self.ax1.plot(self.w, self.extrold + self.tfun[0],"--",color="y")
        self.ax1.plot(self.w, self.objpos + self.tfun[0],"--",color="b")

    def display_redshift(self):
        dw = (self.w1-self.w0)/200.
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

    def bounds(self):
        sigma = 5

        vstd = 1.4826*MAD(self.big.flat)
        vmed = np.median(self.big.flat)
        ystd = 1.4826*MAD(self.spec1d[spec_i,:])
        ymed = np.median(self.spec1d[spec_i,:])

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
            p = fig2.add_subplot(111)
            p.hist(new_z,bins)
            p.text(0.0,1.025,"N = %i" % (len(new_z)), transform=p.transAxes)
            p.text(0.9,1.025,"zq > %i " % (zq_th), transform=p.transAxes)
            p.set_xlabel("Redshift z")
            plt.show()

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
            
        if event.key in ["[","]"] and wlist and self.redshift:
            self.z = self.setup_obj["z"]
            if type(self.z) is NoneType: self.z = 0

            if event.key == "[": self.z -=0.01
            if event.key == "]": self.z +=0.01

            objD[rows[self.fr]]["setup"]["z"] = self.z

        if event.key in ["{","}"] and wlist and self.redshift:
            self.z = self.setup_obj["z"]
            if type(self.z) is NoneType: self.z = 0

            if event.key == "{": self.z -=0.001
            if event.key == "}": self.z +=0.001

            objD[rows[self.fr]]["setup"]["z"] = self.z

        # Redshift Quality
        if event.key == "0":
            self.zq = 0
            objD[rows[self.fr]]["setup"]["z"] = self.zq
        if event.key == "1":
            self.zq = 1
            objD[rows[self.fr]]["setup"]["z"] = self.zq
        if event.key == "2":
            self.zq = 2
            objD[rows[self.fr]]["setup"]["z"] = self.zq
        if event.key == "3":
            self.zq = 3
            objD[rows[self.fr]]["setup"]["z"] = self.zq
        if event.key == "4":
            self.zq = 4
            objD[rows[self.fr]]["setup"]["z"] = self.zq
        if event.key == "5":
            self.zq = 5
            objD[rows[self.fr]]["setup"]["z"] = self.zq

        # Predefined Redshifts
        if event.key == "!":
            self.z = 0.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "@":
            self.z = 0.10
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "#":
            self.z = 0.35
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "$":
            self.z = 0.75
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "%":
            self.z = 1.50
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "^":
            self.z = 2.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "&":
            self.z = 3.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "*":
            self.z = 4.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == "(":
            self.z = 5.00
            objD[rows[self.fr]]["setup"]["z"] = self.z
        if event.key == ")":
            self.z = 6.00
            objD[rows[self.fr]]["setup"]["z"] = self.z


        if event.key == "C":
            try:
                self.comment = raw_input("Enter comment? ")
                objD[rows[self.fr]]["setup"]["comment"] = self.comment 
            except:
                print "Error!"

        if event.key == 'p':

            # cols = ["ra","dec","object","type","z","zq","flag","features", "comment"]

            cols = ["ra","dec","object","z","zq","comment"]
            fmt  = ["%s","%s","%-18s","%7.4f","%6i","%s"]

            for i in objD:
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

try:
    ph1 = os.getenv("PYTHONHOME1")
    if not ph1: raise "PYTHONHOME1 environmental variable not set!"
    wavefile = "/".join([ph1,"/datafiles/linelists/galaxylines.dat"])
except:
    wavefile = None

# raw
spec_i = 0
noise_i = 1
# fluxed
spec_i = 4
noise_i = 5
# comparison
spec_i = 0
noise_i = 4 # fluxed spectra



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

smf = args.smf
fdir = args.dir
fr = args.fr
trace = args.trace




db = simpledb(smf,fdir)
#print db
objD = db.dict["objects"]
rows = db.row_search()
#print rows
#for d in objD: print d,objD[d]

N = len(rows)


### linelists 
if wavefile and os.path.exists(wavefile):
    wlines = open(wavefile,"r").readlines()
    wlines = filter(lambda l: l[0] != "#", wlines)
    wlist = map(str.split, wlines)
else:
    wlist = None


# import traces
tfun = None
ind = 0
if trace:
    tfun = cPickle.load(open(trace))
    # Compatability mode (older version)
    #if type(tfun) is ListType: ind = 0
    # Temporary mode, future mode to be added to fs
    #if type(tfun) is np.ArrayType:
    #    if tfun.shape[0] > 1: ind = 1
    #    elif tfun.shape[0] == 1: ind = 0
    #print type(tfun)
    #print tfun.shape

#objD["trace"] = tfun
#print objD["trace"]
db.dict["trace"] = tfun





# reset matplotlib standard keymap
for i in plt.rcParams:
    if "keymap" in i:
        plt.rcParams[i] = ''

#sys.exit()

fig = plt.figure()
fig.canvas.set_window_title('g_spec.py')
#fig = plt.figure(figsize=(10,8))
SO = SingleObject(fig,fr,N)
SO.connect()
plt.show()

db.write()


