from scipy.ndimage import median_filter
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

def binnedArray(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[(i)*nbin:(i+1)*nbin],axis=0) for i in range(int(len(arr)/nbin)) ])

def binnedArray2(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[:,(i)*nbin:(i+1)*nbin],axis=1) for i in range(int(arr.shape[1]/nbin)) ]).T

def estimate_S(X, L, offposFlag, k=100, s=5, pos_only=True):
    S = X - L
    S[offposFlag] = 0.0
    spec = median_filter(np.sum(S,axis=0), s)

    if pos_only:
        spec = np.minimum(spec, 0)

    spec = np.abs(spec)

    for i in (range(len(spec)-k)):
        ch = np.nanargmin(spec)
        S[:,ch] = 0.0
        spec[ch] = np.nan

    return S


def estimate_L(X, S, n=10):
    R = X - S
    R0 = np.mean(R, axis=0)

    model = TruncatedSVD(n)
    C = model.fit_transform(R - R0)
    P = model.components_
    L = np.full_like(X, C @ P) + R0

    return L

class GoDecDataSet:
    def __init__(self, B4Rdataset, binning=1, pol=['pol1'], sideband='lsb', globLogDir='.'):

        usespw = [ i for i in range(B4Rdataset.nIF) if (B4Rdataset.sbdef[i] == sideband) and ( B4Rdataset.poldef[i] in pol) ]

        self.B4Rdataset = B4Rdataset
        self.binning = binning
        self.T_amb = 273.+B4Rdataset.Temperature_cal
        self.usespw = usespw
        self.pol = pol
        self.sideband = sideband
        self.lsbfreq  = binnedArray(B4Rdataset.flsb, binning)
        self.usbfreq  = binnedArray(B4Rdataset.fusb, binning)
        self.obsnum   = B4Rdataset.obsnum
        self.globLogDir = globLogDir

    def B4Rdata2GoDecInput(self):
        rspec = self.B4Rdataset.rspec

        self.offposFlag = np.array([x == ('REF', 'S') for x in self.B4Rdataset.xffts.scanPattern])

        ndat = self.B4Rdataset.xffts.getSize()
        dP = []
        for i in tqdm(np.arange(ndat)) :
            Tint = self.B4Rdataset.xffts.binHead(i)['integTime']
            spec = binnedArray2(self.B4Rdataset.xffts.readData(i)/Tint-rspec, self.binning)
            dSignal = np.mean(spec[self.usespw,:], axis=0)
            dP.append(dSignal)

        self.ndat = ndat
        self.dP = np.array(dP)
        self.X  = np.log(-self.dP/self.T_amb)


    def GoDecCal(self,k=25,s=5):

        S = np.zeros_like(self.X)

        print('### GoDec Calibration with k='+str(k)+', s='+str(s))
        for i in tqdm(range(30)):
            L = estimate_L(self.X, S, n=5)
            S = estimate_S(self.X, L, self.offposFlag, k=k, s=s, pos_only=True)

        T_cal = self.T_amb * (1-np.exp(self.X-L))
        T_cal = T_cal[~self.offposFlag]

        self.T_cal = T_cal

    def outputResult(self,plot=True):

        if self.sideband == 'lsb':
            freq = self.lsbfreq
        elif self.sideband == 'usb':
            freq = self.usbfreq

        spec = np.mean(self.T_cal,axis=0)
        os.system('mkdir -p '+self.globLogDir+'/'+str(self.obsnum))
        os.system('mkdir -p '+self.globLogDir+'/'+str(self.obsnum)+'/GoDec')
        pols = ''
        for i in self.pol:
            pols = pols+i
        basename = self.globLogDir+'/'+str(self.obsnum)+'/GoDec/'+str(self.obsnum).zfill(6)+'.GoDec.binning_'+str(self.binning)+'.sideband_'+self.sideband+'.pol_'+pols
        np.savez(basename+'.npz',freq=freq,spec=spec)

        if plot:
            plt.close()
            plt.step(freq,spec,'b',where='mid')
            plt.xlabel('frequency [GHz]')
            plt.ylabel(r'T$_\mathrm{a}^*$ [K] (GoDec)')
            plt.savefig(basename+'.png')
            plt.close()
