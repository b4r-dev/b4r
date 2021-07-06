#!/usr/bin/env python3
import b4rpipe.LibB4Rtools as Lib
import b4rpipe.b4rquery as b4rQ
import b4rpipe.godec as godec
import b4rpipe as Bp
#import importlib
import os
import numpy as np
#importlib.reload(Lib)

Bp.globBaseDir = './rawdata'
Bp.globLogDir = './calibrated'

def PipelineAnalysis(obsnum,dAZ=0.,dEL=0.,dAZdELmode=False,DataDownload=False,username='',password=''):

    Lib.globBaseDir = Bp.globBaseDir
    Lib.globLogDir = Bp.globLogDir

    os.system('mkdir -p '+Lib.globLogDir)
    os.system('mkdir -p '+Lib.globLogDir+'/'+str(obsnum))

    logf = open(Lib.globLogDir+'/'+str(obsnum)+'/PipelineAnalysis.'+str(obsnum)+'.log','w')

    if DataDownload:
        os.system('mkdir -p '+Bp.globBaseDir)
        Query = b4rQ.B4Rquery(username=username,password=password,localBaseDir=Bp.globBaseDir)
        Query.SearchAndDownload(obsnum)
        Query.SearchAndDownload(obsnum-1)

    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=8)
        obj.LinePointingQlook()
        logf.write('LinePointing: PASS'+'\n')
    except:
        logf.write('LinePointing: FAIL'+'\n')

    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=512,noRefCal=True)
        obj.ContPointingQlook()
        logf.write('ContPointing: PASS'+'\n')
    except:
        logf.write('ContPointing: FAIL'+'\n')

    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=8)
        obj.PsQlook(highz=False)
        logf.write('Ps(nearby): PASS'+'\n')
    except:
        logf.write('Ps(nearby): FAIL'+'\n')

    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=256)
        obj.PsQlook(highz=True)
        logf.write('Ps(High-z): PASS'+'\n')
    except:
        logf.write('Ps(High-z): FAIL'+'\n')


    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=1)
        for pol in ['pol1','pol2']:
            for sideband in ['usb','lsb']:
                objGoDec = godec.GoDecDataSet(obj, binning=256, pol=[pol], sideband=sideband, globLogDir=Lib.globLogDir)
                objGoDec.B4Rdata2GoDecInput()
                objGoDec.GoDecCal()
                objGoDec.outputResult()
        logf.write('GoDec(High-z): PASS'+'\n')
    except:
        logf.write('GoDec(High-z): FAIL'+'\n')
    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=1)
        obj.createPsData()
        logf.write('PsData: PASS'+'\n')
    except:
        logf.write('PsData: FAIL'+'\n')

    try:
        obj = Lib.B4Rdataset(obsnum=obsnum,calnum=obsnum-1)
        obj.Pipeline(binning=1,noRefCal=True)
        obj.createMS2(dAZ=dAZ,dEL=dEL,dAZdELmode=dAZdELmode)
        logf.write('MS2: PASS'+'\n')
    except:
        logf.write('MS2: FAIL'+'\n')

    logf.close()

def PipelineAnalysisBatchRun():
    print('### Batch run start!! ###')
    obsnum_list = (np.array(Lib.returnFileList(mode='OBS'))[:,0]).astype('int')
    obsnum_list = obsnum_list[obsnum_list>79999] # 2019

    for obsnum in obsnum_list:
        PipelineAnalysis(obsnum)

def debug():
    print(Bp.globBaseDir)
    print(Bp.globLogDir)
