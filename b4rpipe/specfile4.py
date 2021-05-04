import numpy as np
import struct
import os
import glob
import math

from socket import *
import time
#from thread import *
import os
import subprocess

#from pyqtgraph.Qt import QtGui, QtCore
#import pyqtgraph as pg
#import pyqtgraph.exporters




class SpecFile :
    maxch = 32768
    headSize = 256
    unitSize = headSize + maxch*4

    def __init__ (self, name, directory = '', spw=[0,1,2,3], chrange=[0,None]) :
        self.name = name
        self.directory = directory
        self.fullPath = '%s%s'%(directory, name)
        self.narray  = len(spw)
        self.fnames = []
        self.scanPattern  = []
        self.scanIds      = []
        self.spw  = spw

        if chrange[0] == None : chrange[0] = SpecFile.maxch
        if chrange[1] == None : chrange[1] = SpecFile.maxch
        if chrange[0] < 0 : chrange[0] = chrange[0]+SpecFile.maxch
        if chrange[1] < 0 : chrange[1] = chrange[1]+SpecFile.maxch
        self.chrange = sorted(chrange)
        self.nchan  = self.chrange[1] - self.chrange[0]

    def getSize(self) :
        size_raw = os.stat(self.fullPath + '.01').st_size
        if size_raw < SpecFile.unitSize : return 0
        return int(math.floor(size_raw / SpecFile.unitSize))

    def binHead(self, idx) :
        fp = open(self.fullPath+'.01', 'rb')
        fp.seek(idx * SpecFile.unitSize)
        bindata = fp.read(SpecFile.headSize)
        bin_date   = bindata[:28]
        bin_obsNum = bindata[32:40]
        bin_bufPos = bindata[40:48]
        bin_scanMode =bindata[48:56]
        bin_chopper  =bindata[56:64]
        bin_scanCount=bindata[64:72]
        bin_specCount=bindata[72:80]
        bin_integTime=bindata[80:88]
        return { 'date' : bin_date,
                 'obsNum'    : struct.unpack_from('<l', bin_obsNum)[0],
                 'scanCount' : struct.unpack_from('<l', bin_scanCount)[0],
                 'specCount' : struct.unpack_from('<l', bin_specCount)[0],
                 'bufPos'    : bin_bufPos  [:bin_bufPos .find(b'\x00')].decode(),
                 'chopperPos': bin_chopper [:bin_chopper.find(b'\x00')].decode(),
                 'scanMode'  : bin_scanMode[:bin_scanMode.find(b'\x00')].decode(),
                 'integTime' : struct.unpack_from('<l', bin_integTime)[0]};

    def analyScanPattern(self) :
        self.scanPattern = []
        self.chopPattern = []
        self.scanIds     = []
        self.calCount    = []
        self.onCount     = []
        self.offCount    = []

        theSize = self.getSize()
        i_cal, i_on, i_off = -1,-1,-1
        for idx in range(theSize) :
            theHead = self.binHead(idx)
            self.scanPattern.append((theHead['bufPos'], theHead['chopperPos']))
            self.scanIds    .append(theHead['scanCount'])

            isnewScan = ( self.scanIds[-1] != self.scanIds[-2] ) if idx > 1 else True
            thisScan = self.scanPattern[-1]

            if thisScan[0] == 'CAL' and isnewScan : i_cal += 1
            if thisScan[0] == 'ON'  and isnewScan : i_on  += 1
            if thisScan[0] == 'REF' and isnewScan : i_off += 1

            self.calCount.append(i_cal)
            self.onCount .append(i_on)
            self.offCount.append(i_off)

        self.numCal = max(self.calCount)+1 if len(self.calCount)> 0  else 0
        self.numOn  = max(self.onCount) +1 if len(self.onCount)> 0  else 0
        self.numOff = max(self.offCount)+1 if len(self.offCount)> 0  else 0
        self.numScan= max(self.scanIds)+1 if len(self.scanIds)> 0  else 0

    def readSpec(self, idx):

        data = [ None ] * (self.narray)
        date = ''

#        for i in range(self.narray) :
        for i in range(self.narray) :
            fp = open(self.fullPath+'.%02d'%(self.spw[i]+1), 'rb')
            fp.tell()
            fp.seek(idx * SpecFile.unitSize)
            bindata = fp.read(SpecFile.headSize )#unitSize)

            bin_date   = bindata[:28]
            bin_obsNum = bindata[32:40]
            bin_bufPos = bindata[40:48]
            bin_scanMode =bindata[48:56]
            bin_chopper  =bindata[56:64]
            bin_scanCount=bindata[64:72]
            bin_specCount=bindata[72:80]
            bin_integTime=bindata[80:88]

#            data[i] =  [ struct.unpack_from('<f', bindata, 256 + 4*k)[0] for k in range(32768)]
#            data[i] =  [ struct.unpack_from('<f', bindata, 256 + 4*k)[0] for k in range(self.chrange[0], self.chrange[1])]
#            data[i] = np.frombuffer(bindata, offset = 256 + self.chrange[0]*4, count=self.nchan, dtype='f') #dtype=float)

            fp.seek(idx * SpecFile.unitSize + 256 + self.chrange[0]*4)
            data[i]  = np.fromfile(fp, count=self.nchan, dtype='f') #dtype=float)

        return { 'date' : bin_date,
                 'data' : np.array(data),
                 'obsNum'    : struct.unpack_from('<l', bin_obsNum)[0],
                 'scanCount' : struct.unpack_from('<l', bin_scanCount)[0],
                 'specCount' : struct.unpack_from('<l', bin_specCount)[0],
                 'bufPos'    : bin_bufPos  [:bin_bufPos .find(b'\x00')].decode(),
                 'chopperPos': bin_chopper [:bin_chopper.find(b'\x00')].decode(),
                 'scanMode'  : bin_scanMode[:bin_scanMode.find(b'\x00')].decode(),
                 'integTime' : struct.unpack_from('<l', bin_integTime)[0]        };

    def readSpec0(self, idx):

        data = [ None ] * (self.narray)
        date = ''

#        for i in range(self.narray) :
        for i in range(self.narray) :
            fp = open(self.fullPath+'.%02d'%(self.spw[i]+1), 'rb')
            fp.tell()
            fp.seek(idx * SpecFile.unitSize)
            bindata = fp.read(SpecFile.unitSize)

            bin_date   = bindata[:28]
            bin_obsNum = bindata[32:40]
            bin_bufPos = bindata[40:48]
            bin_scanMode =bindata[48:56]
            bin_chopper  =bindata[56:64]
            bin_scanCount=bindata[64:72]
            bin_specCount=bindata[72:80]
            bin_integTime=bindata[80:88]

#            data[i] =  [ struct.unpack_from('<f', bindata, 256 + 4*k)[0] for k in range(32768)]
#            data[i] =  [ struct.unpack_from('<f', bindata, 256 + 4*k)[0] for k in range(self.chrange[0], self.chrange[1])]
            data[i] = np.frombuffer(bindata, offset = 256 + self.chrange[0]*4, count=self.nchan, dtype='f') #dtype=float)

#            fp.seek(idx * SpecFile.unitSize + 256 + self.chrange[0]*4)
#            data[i]  = np.fromfile(fp, count=self.nchan, dtype='f') #dtype=float)

        return { 'date' : bin_date,
                 'data' : np.array(data),
                 'obsNum'    : struct.unpack_from('<l', bin_obsNum)[0],
                 'scanCount' : struct.unpack_from('<l', bin_scanCount)[0],
                 'specCount' : struct.unpack_from('<l', bin_specCount)[0],
                 'bufPos'    : bin_bufPos  [:bin_bufPos .find(b'\x00')].decode(),
                 'chopperPos': bin_chopper [:bin_chopper.find(b'\x00')].decode(),
                 'scanMode'  : bin_scanMode[:bin_scanMode.find(b'\x00')].decode(),
                 'integTime' : struct.unpack_from('<l', bin_integTime)[0]        };




    def readData(self, idx):

        data = [ None ] * (self.narray)

        for i in range(self.narray) :
            fp = open(self.fullPath+'.%02d'%(self.spw[i]+1), 'rb')
            fp.tell()
            fp.seek(idx * SpecFile.unitSize + 256 + self.chrange[0]*4)
            data[i]  = np.fromfile(fp, count=self.nchan, dtype='f') #dtype=float)

        return np.array(data)


    def getLastAveragedSpec(self, scanType, chopPos) :
        self.analyScanPattern()
        theNum = len(self.scanPattern)
        scanPattern_r = copy(self.scanPattern)[::-1].tolist()

        if not ( [scanType, chopPos] in scanPattern_r ) :
            return None

        theLastIdx = scanPattern_r.index([scanType, chopPos])
        spec2avg = []
        print(theLastIdx)
        for i in range(theLastIdx, theNum) :
            if scanPattern_r[i] != [scanType, chopPos] : break
            print(i, self.scanPattern[theNum-1-i])
            sp = self.readSpec( theNum -1 - i )
            spec2avg.append(sp['data'])

        print(len(spec2avg))
        return np.average(spec2avg, axis= 0)

    def getAveragedSpec(self, scanId) :
        self.analyScanPattern()
        theNum = len(self.scanIds)
        spec2avg = np.array([ self.readSpec(i)['data'] for i in range(theNum) if self.scanIds[i] == scanId ])
#        intgTime = ones(spec2avg.shape[:2])
        intgTime = np.array([self.binHead(i)['integTime']*1.e-6 for i in range(theNum) if self.scanIds[i] == scanId ])
        print(intgTime.shape)
#        for i in range(len(intgTime)) :
#            spec2avg[i] *= intgTime[i]
        sumT = np.sum(intgTime)
        print(sumT)
        return np.sum(spec2avg, axis=0)/ sumT
#        return array( [ spec2avg[i]*intgTime[i] for i in range(len(spec2avg))] )

    def getAveragedSpecDiv(self, scanId, n, N) :
        self.analyScanPattern()
        theNum = len(self.scanIds)
        spec2avg = [ self.readSpec(i)['data'] for i in range(theNum) if self.scanIds[i] == scanId ]
        print( [ (self.scanIds[i], self.scanPattern[i]) for i in range(theNum) if self.scanIds[i] == scanId ])
        print(len(spec2avg))
        ntot= len(spec2avg)
        nunit = ntot/N
        st = n * nunit
        ed = st + nunit + 1
        return np.average(spec2avg[st:ed], axis= 0)

    def getLastOnOff(self) :
        lastOn  = self.getLastAveragedSpec('ON' , 'S')
        lastOff = self.getLastAveragedSpec('REF', 'S')
        if (lastOn == None) or (lastOff == None) : return None
        else : return lastOn - lastOff

    def getLastRsky(self) :
        lastR = self.getLastAveragedSpec('CAL', 'R')
        lastS = self.getLastAveragedSpec('CAL', 'S')
        if (lastR == None) or (lastS == None) : return None
        else : return lastR - lastS

    def getLatCal(self) :
        return 300 * self.getLastOnOff() / self.getLastRsky()

    def getOnOff(self, onIdx) :
        self.analyScanPattern()
        offIdx = self.offCount[ self.onCount.index(onIdx) ]
        if offIdx < 0 : return None

        specOn  = self.getAveragedSpec(self.scanIds[self.onCount.index(onIdx)])
        specOff = self.getAveragedSpec(self.scanIds[self.offCount.index(offIdx)])

        return specOn - specOff

    def getOnOffDiv(self, onIdx, n, N ) :
        self.analyScanPattern()
        offIdx = self.offCount[ self.onCount.index(onIdx) ]
        if offIdx < 0 : return None

        specOn  = self.getAveragedSpecDiv(self.scanIds[self.onCount.index(onIdx)], n, N)
        specOff = self.getAveragedSpec(self.scanIds[self.offCount.index(offIdx)])

        return specOn - specOff

    def getSingleOnOff(self, onIntegIdx) :
        self.analyScanPattern()
        idx = 0
        onCount = 0
        isFound = False
        for x in self.scanPattern :
            if x[0] == 'ON' :
                if onCount == onIntegIdx :
                    isFound = True
                    break
                else : onCount += 1
            idx +=1

        if not isFound : return None

        theSpec = self.readSpec(idx)
        offIdx = self.offCount[ idx ]
        onIdx = self.onCount[ idx ]
        if offIdx < 0 : return None

        print(idx, offIdx, onIdx)
        specOn = theSpec['data'] / (float(theSpec['integTime'])*1.e-6)
        specOff= self.getAveragedSpec(self.scanIds[self.offCount.index(offIdx)])

        return (specOn - specOff), [idx,offIdx,onIdx]
