
import glob,re,os,shutil,sys
import cPickle

import numpy as np

class simpledb:
    def __init__(self,smf,fdir):

        self.smf = smf
        self.fdir = fdir

        self.mask = self.smf.replace(".SMF","")

        L_all = open(self.smf,"r").readlines()
        L = [l.split() for l in L_all]

        self.smf_type = [l[0] for l in L if l[0] in ["HOLE","SLIT"]]
        self.smf_ids  = [l[1] for l in L if l[0] in ["HOLE","SLIT"]]
        self.ra       = [l[2] for l in L if l[0] in ["HOLE","SLIT"]]
        self.dec      = [l[3] for l in L if l[0] in ["HOLE","SLIT"]]
        self.smf_x    = np.array(map(float,[l[4] for l in L if l[0] in ["HOLE","SLIT"]]))
        self.smf_y    = np.array(map(float,[l[5] for l in L if l[0] in ["HOLE","SLIT"]]))

        self.refobjs  = [l[1] for l in L if l[0] == "HOLE"]

        self.open()
        if not self.dict.has_key("objects") or not self.dict.has_key("smf"): self.create()
        else: self.update()

    def open(self):
        if os.path.exists(self.mask+".p"):
            # if file is not 0 bytes
            if os.path.getsize(self.mask+".p"):
                shutil.copyfile(self.mask+".p", self.mask+".p.bkup")
            else:
                #EOFError
                print "Python pickle file is corrupted! Use pickle backup!"
                print "Exiting..."
                sys.exit()
            self.dict = cPickle.load(open( self.mask+".p", "rb" ))
        else:
            self.dict = {}

    def write(self):
        cPickle.dump(self.dict, open(self.mask+".p", "wb" ))

    def create(self):
        print "no db"
        self.dict["objects"] = {}
        self.dict["smf"] = {"type":self.smf_type,
                            "ids":self.smf_ids,
                            "x":self.smf_x,
                            "y":self.smf_y}

        r2 = self.dict["objects"]
        for i,j in enumerate(self.smf_ids):

            r2[i] = {}
            r2[i]["setup"] = {"object":j,
                              "smf_type":self.smf_type[i],
                              "smf_x":self.smf_x[i],
                              "smf_y":self.smf_y[i],
                              "ra":self.ra[i],
                              "dec":self.dec[i],
                              #"spec2d":self.fdir+"/"+j+"sum.fits",
                              #"spec2d_sig":self.fdir+"/"+j+"sumsig.fits",
                              #"spec2d_sky":self.fdir+"/"+j+"sumsky.fits",
                              #"spec2d_ext":self.fdir+"/"+j+"sumext.fits",
                              #"spec1d":self.fdir+"/"+j+"_1dspec.fits",
                              "spectrum":"",
                              "type":"",
                              "features":"",
                              "z":None,
                              "zq":None,
                              "flag":None,
                              "comment":""}

            spec_key = ["spec2d","spec2d_sig","spec2d_sky","spec2d_ext"]
            spectra = ["sum.fits","sumsig.fits","sumsky.fits","sumext.fits"]
            spec1d = self.fdir+"/"+j+"_1dspec.fits"
            spec1d_flux = self.fdir+"/"+j+"_1dspec_flux.fits"

            for k,s in zip(spec_key,spectra):
                f = self.fdir+"/"+j+s
                if not os.path.exists(f):
                    print "Missing file %s!!!" % f
                    f = None
                    #spec1d = None

                r2[i]["setup"][k] = f

            if not os.path.exists(spec1d):
                spec1d = None
            if not os.path.exists(spec1d_flux):
                spec1d_flux = None

            r2[i]["setup"]["spec1d"] = spec1d
            r2[i]["setup"]["spec1d_flux"] = spec1d_flux


    def update(self):
        r2 = self.dict["objects"]
        for i,j in enumerate(self.smf_ids):

            SetupDict = {"object":j,
                         "smf_type":self.smf_type[i],
                         "smf_x":self.smf_x[i],
                         "smf_y":self.smf_y[i],
                         "ra":self.ra[i],
                         "dec":self.dec[i],
                         #"spec2d":self.fdir+"/"+j+"sum.fits",
                         #"spec2d_sig":self.fdir+"/"+j+"sumsig.fits",
                         #"spec2d_sky":self.fdir+"/"+j+"sumsky.fits",
                         #"spec2d_ext":self.fdir+"/"+j+"sumext.fits",
                         #"spec1d":self.fdir+"/"+j+"_1dspec.fits",
                         "spectrum":"",
                         "type":"",
                         "features":"",
                         "z":None,
                         "zq":None,
                         "flag":None,
                         "comment":""}

            spec_key = ["spec2d","spec2d_sig","spec2d_sky","spec2d_ext"]
            spectra = ["sum.fits","sumsig.fits","sumsky.fits","sumext.fits"]
            spec1d = self.fdir+"/"+j+"_1dspec.fits"
            spec1d_flux = self.fdir+"/"+j+"_1dspec_flux.fits"

            for k,s in zip(spec_key,spectra):
                f = self.fdir+"/"+j+s
                if not os.path.exists(f):
                    print "Missing file %s!!!" % f
                    f = None
                    #spec1d = None

                SetupDict[k] = f

            if not os.path.exists(spec1d):
                spec1d = None
            if not os.path.exists(spec1d_flux):
                spec1d_flux = None

            SetupDict["spec1d"] = spec1d
            SetupDict["spec1d_flux"] = spec1d_flux

            for sd in SetupDict:
                if not r2[i]["setup"].has_key(sd):
                    r2[i]["setup"][sd] = SetupDict[sd]



    def row_search(self):
        return self.dict["objects"].keys()
