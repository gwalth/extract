
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

    def open(self):
        if os.path.exists(self.mask+".p"):
            shutil.copyfile(self.mask+".p", self.mask+".p.bkup")
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
                              "ra":self.ra[i],
                              "dec":self.dec[i],
                              "spec2d":self.fdir+"/"+j+"sum.fits",
                              "spec2d_sig":self.fdir+"/"+j+"sumsig.fits",
                              "spec2d_sky":self.fdir+"/"+j+"sumsky.fits",
                              "spec2d_ext":self.fdir+"/"+j+"sumext.fits",
                              "spec1d":self.fdir+"/"+j+"_1dspec.fits",
                              "type":"",
                              "features":"",
                              "z":None,
                              "zq":None,
                              "flag":None,
                              "comment":""}


    def row_search(self):
        return self.dict["objects"].keys()
