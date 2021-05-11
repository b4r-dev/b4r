
import matplotlib
#matplotlib.use('TkAgg')
import os

import netCDF4 as nc

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime as dt
import astropy.stats
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, LSR, AltAz, ICRS
import astropy.units as u

import importlib
import scipy.optimize
from tqdm import tqdm

from copy import copy
from math import pi
import re
import warnings
warnings.simplefilter('ignore')


from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.auto_max_age = None

# submodules
import b4rpipe.makedatalist as md
import b4rpipe.specfile4 as sp
import b4rpipe.CreateMs2 as ms2

globBaseDir = './rawdata'
globLogDir = './calibrated'

def binnedArray(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[(i)*nbin:(i+1)*nbin],axis=0) for i in range(int(len(arr)/nbin)) ])

def binnedArray2(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[:,(i)*nbin:(i+1)*nbin],axis=1) for i in range(int(arr.shape[1]/nbin)) ]).T

def modelBeam(xymesh, posx, posy, widthx, widthy, peak, sky) : #modelPar) :
    gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
    axis_x = xymesh[0][0]
    axis_y = xymesh[1][:,0]
    modelx = gsbeam( axis_x - posx, widthx)[:,np.newaxis]
    modely = gsbeam( axis_y - posy, widthy)[:,np.newaxis]
    return (np.dot(modelx, modely.T) * peak + sky).flatten()

def checkspec(filtered_data,startch=0,merginch=3):

    ip=0
    while((filtered_data[startch+ip] == 1 or filtered_data[startch+ip+1] == 1) and startch+ip+1<filtered_data.shape[0]):
        ip=ip+1

    im=0
    while((filtered_data[startch-im] == 1 or filtered_data[startch-im-1] == 1) and startch-im-1>=0):
        im=im+1

    return [startch-im-merginch,startch+ip+merginch]

def RMSexpect(Tsys_value,df,tint=10,npol=1):
    return 2**0.5*Tsys_value/(df*1.0e9*tint*npol)**0.5

def AZEL2RADEC(az,
               el,
               utctime,
               lon = -97.31461167440136,
               lat = 18.98607157170612,
               height = 4640.):

	# modules
	from astropy import coordinates
	from astropy.time import Time
	import astropy.units as u

	site = coordinates.EarthLocation.from_geodetic(lon=lon,lat=lat,height=height)
	c = str(az) + ' ' + str(el)

	radec = coordinates.SkyCoord(c, unit=(u.deg,u.deg),frame='altaz',obstime=utctime,location=site).transform_to('icrs')
	return radec.ra.value, radec.dec.value


def returnFileList(date='', source='', mode='', freq='', obsgoal='') :
    utc = dt.timezone(dt.timedelta(hours=+0), 'utc')
    info_list = []

    for x in md.getFileList(globBaseDir, ffilter=date) :
        p = md.getNCprops(x)
        if p == None : continue
        if re.compile(source).search(p['srcname']) == None : continue
        if re.compile(mode)  .search(p['obsmode']) == None : continue
        if re.compile(freq)  .search(str(p['obsfreq'])) == None : continue
        if re.compile(obsgoal)  .search(str(p['obsgoal'])) == None : continue

        info = [[str(p['obsnum']).zfill(6), str(dt.datetime.fromtimestamp(p['tm'][0], tz=utc)), '{:3.1f}'.format(p['tm'][1]-p['tm'][0]), p['srcname'], str(p['obsfreq']), p['obsmode'], p['obsgoal']]]
        info_list = info_list + info

    return info_list

class B4Rdataset :

    def __init__(self,obsnum=-1,calnum=-1):
        self.obsnum = obsnum
        self.calnum = calnum
        print('### Info ')

        if obsnum == -1:
            self.ncfname   = None
            self.befname   = None
            print('obsnum=None')
            print('')

        else:
            print('obsnum='+str(obsnum))
            self.ncfname = md.getLMTtpmByObsNum(obsnum, globBaseDir)
            self.befname = md.getXFFTSbyObsNum(obsnum, globBaseDir)

        if calnum == -1:
            self.calncfname = None
            self.calbename = None
            print('calnum=None')
            print('')
        else:
            print('calnum='+str(calnum))
            self.calncfname =  md.getLMTtpmByObsNum(calnum, globBaseDir)
            self.calbename = md.getXFFTSbyObsNum(calnum, globBaseDir)


    def loadlmttpm(self) :

        if self.ncfname != None:

            ncfile= nc.Dataset(self.ncfname)
            ncprop = md.getNCprops(self.ncfname)

            utc = dt.timezone(dt.timedelta(hours=+0), 'utc')

            nctime = ncfile.variables['Data.TelescopeBackend.TelTime'][:]
            bufpos = ncfile.variables['Data.TelescopeBackend.BufPos'][:]

            self.az_map = ncfile['Data.TelescopeBackend.TelAzMap'][:]
            self.el_map = ncfile['Data.TelescopeBackend.TelElMap'][:]
            self.az_sky = ncfile['Data.TelescopeBackend.TelAzSky'][:]
            self.el_sky = ncfile['Data.TelescopeBackend.TelElSky'][:]
            self.source_az = ncfile['Data.TelescopeBackend.SourceAz'][:]
            self.source_el = ncfile['Data.TelescopeBackend.SourceEl'][:]

            self.fif    = [ncfile['Header.B4r.If1Freq'][:][0], ncfile['Header.B4r.If2Freq'][:][0]]
            self.fline  = ncfile['Header.B4r.LineFreq'][0]

            self.srcname = ncprop['srcname']
            self.source_ra = ncfile['Header.Source.Ra'][0]
            self.source_dec = ncfile['Header.Source.Dec'][0]
            self.Pid = ncprop['projid']
            self.Observer =  ncprop['observer']
            self.sysvel = ncfile['Header.Source.Velocity'][0]

            self.Temperature = ncfile['Header.Weather.Temperature'][0]
            self.tau = ncfile['Header.Radiometer.Tau'][0]
            self.WindSpeed1 = ncfile['Header.Weather.WindSpeed1'][0]
            self.WindDir1 = ncfile['Header.Weather.WindDir1'][0]
            self.WindSpeed2 = ncfile['Header.Weather.WindSpeed2'][0]
            self.WindDir2 = ncfile['Header.Weather.WindDir2'][0]

            self.dAZuser = ncfile['Header.PointModel.AzUserOff'][0]
            self.dELuser = ncfile['Header.PointModel.ElUserOff'][0]
            self.dAZModel = ncfile['Header.PointModel.AzPointModelCor'][0]
            self.dELModel = ncfile['Header.PointModel.ElPointModelCor'][0]
            self.dAZM2Cor = ncfile['Header.PointModel.AzM2Cor'][0]
            self.dELM2Cor = ncfile['Header.PointModel.ElM2Cor'][0]


            self.antonpos = bufpos == 0
            self.antoffpos= bufpos == 1

            self.t_ant = np.array([dt.datetime.fromtimestamp(t, tz=utc) for t in nctime], 'datetime64[us]')

            ncfile.close()

        if self.calncfname != None:

            ncfile= nc.Dataset(self.calncfname)
            ncprop = md.getNCprops(self.calncfname)

            utc = dt.timezone(dt.timedelta(hours=+0), 'utc')

            nctime = ncfile.variables['Data.TelescopeBackend.TelTime'][:]
            bufpos = ncfile.variables['Data.TelescopeBackend.BufPos'][:]

            self.fif_cal    = [ncfile['Header.B4r.If1Freq'][:][0], ncfile['Header.B4r.If2Freq'][:][0]]
            self.fline_cal  = ncfile['Header.B4r.LineFreq'][0]
            self.srcname_cal = ncprop['srcname']

            self.Temperature_cal = ncfile['Header.Weather.Temperature'][0]
            self.tau_cal = ncfile['Header.Radiometer.Tau'][0]
            self.WindSpeed1_cal = ncfile['Header.Weather.WindSpeed1'][0]
            self.WindDir1_cal = ncfile['Header.Weather.WindDir1'][0]
            self.WindSpeed2_cal = ncfile['Header.Weather.WindSpeed2'][0]
            self.WindDir2_cal = ncfile['Header.Weather.WindDir2'][0]

            self.antonpos_cal = bufpos == 0
            self.antoffpos_cal= bufpos == 1

            self.t_ant_cal = np.array([dt.datetime.fromtimestamp(t, tz=utc) for t in nctime], 'datetime64[us]')

            ncfile.close()


    def loadxffts(self,
                  nIF = 4,
                  sbdef  = [ 'lsb', 'usb', 'lsb', 'usb' ],
                  poldef = [ 'pol1', 'pol1', 'pol2', 'pol2' ],
                  useSB = ['lsb','usb'],
                  usePol= ['pol1','pol2'],
                  timediffTorr= 1) :

        ##############
        # read xffts
        print('### loading xffts data ... ')
        if self.befname != None :

            fif    = self.fif
            fline = self.fline

            df_org = 2.5/32768        #GHz
            iffreq_ch_org = fif[1] + np.arange(0,-2.5, -df_org)
            fusb = fline + iffreq_ch_org
            flsb = fline - iffreq_ch_org

            if 'lsb' in useSB:
                print( "   [OBS] freq range (lsb) : {f} GHz".format(f=[flsb[ 0], flsb[-1]]))
            if 'usb' in useSB:
                print( "   [OBS] freq range (usb) : {f} GHz".format(f=[fusb[-1], fusb[0]]))

            usespw = [ i for i in range(nIF) if (sbdef[i] in useSB ) and ( poldef[i] in usePol) ]

            xfftsf = sp.SpecFile(self.befname, spw=usespw)
            xfftsf.analyScanPattern()

        if self.calbename != None :

            fif_cal    = self.fif_cal
            fline_cal = self.fline_cal

            df_org = 2.5/32768        #GHz
            iffreq_ch_org_cal = fif_cal[1] + np.arange(0,-2.5, -df_org)
            fusb_cal = fline_cal + iffreq_ch_org_cal
            flsb_cal = fline_cal - iffreq_ch_org_cal

            if 'lsb' in useSB:
                print( "   [CAL] freq range (lsb) : {f} GHz".format(f=[flsb_cal[ 0], flsb_cal[-1]]))
            if 'usb' in useSB:
                print( "   [CAL] freq range (usb) : {f} GHz".format(f=[fusb_cal[-1], fusb_cal[0]]))

            usespw_cal = [ i for i in range(nIF) if (sbdef[i] in useSB ) and ( poldef[i] in usePol) ]

            xffts_cal = sp.SpecFile(self.calbename, spw=usespw_cal)
            xffts_cal.analyScanPattern()
            ndat_cal = xffts_cal.getSize()
            betime_cal  = np.array([ xffts_cal.binHead(i)['date'][:-4] for i in range(ndat_cal) ],'datetime64[us]')


        ################################
        # load header inf
        ################################

        if self.befname != None :

            ndat = xfftsf.getSize()
            betime     = np.array([ xfftsf.binHead(i)['date'][:-4] for i in range(ndat) ],'datetime64[us]')   #+  np.timedelta64(9, 'h')
            integtime  = np.array([ xfftsf.binHead(i)['integTime'] for i in range(ndat) ],'timedelta64[us]')
            bebuf  = np.array([ xfftsf.binHead(i)['bufPos'] for i in range(ndat) ])
            beonpos  = bebuf == 'ON'
            beoffpos = bebuf == 'REF'

            ntpsync = np.array( [ xfftsf.binHead(i)['date'][-4:-1] == b'GPS' for i in range(ndat) ])

            tm_PCsync = copy(betime)
            tm_PCsync[ ntpsync == False ] -=  np.timedelta64(9, 'h')

            blnktime = np.timedelta64(1000,'us')

            ntpsync_work        = copy(ntpsync)

            print('### BE time corretion')
            i_needscorrect, = np.where(ntpsync_work==False)
            print('  step1')
            for i in tqdm(i_needscorrect) :
                if i<1 : continue    #just in case
                tmincr = tm_PCsync[i-1] + integtime[i-1] + blnktime
                tmdiff = abs((tmincr-tm_PCsync[i]).astype(float)) * 1.e-6    #microsec->sec
                if tmdiff < timediffTorr and ntpsync_work[i-1] == True :
                    tm_PCsync[i] = tmincr
                    ntpsync_work[i] = True

            i_needscorrect, = np.where(ntpsync_work==False)
            #try:
            print('  step2')
            for i in tqdm(i_needscorrect[::-1]) :
                if i> ndat -1 : continue
                tmincr = tm_PCsync[i+1] - integtime[i] - blnktime
                tmdiff = abs((tmincr-tm_PCsync[i]).astype(float))*1.e-6
                if tmdiff< timediffTorr and ntpsync_work[i+1] == True :
                    tm_PCsync[i]= tmincr
                    ntpsync_work[i] = True
            print('  ')


            ncorrupted = ntpsync.tolist().count(False)
            nrecovered = ((ntpsync == False) & (ntpsync_work == True)).tolist().count(True)
            print(f'   num of corrupted timestamps = {ncorrupted}/{ndat}')
            print(f'   recovered timestamps        = {nrecovered}/{ncorrupted}')
            print( '   flagged time                = {f}%'.format( f = np.sum( integtime[ntpsync_work==False])/np.sum(integtime)))

        if self.befname != None :
            self.nbedat = ndat
            self.betime = tm_PCsync
            self.integtime = integtime
            self.beonpos   = beonpos  & ntpsync_work
            self.beoffpos  = beoffpos & ntpsync_work
            self.corrupted = ntpsync_work==False
            self.xffts      = xfftsf
            self.flsb      = flsb
            self.fusb      = fusb
            #self.f_use      = useFreq
            self.p_use      = usePol
            self.ntpsync = ntpsync
            self.ntpsync_work = ntpsync_work


        if self.calbename != None :
            self.xffts_cal  = xffts_cal
            self.betime_cal = betime_cal

        self.sbdef = sbdef
        self.poldef = poldef
        self.nIF = nIF

    def TsysCal(self,binning=1):
        print(f'   using calib BE data {self.calbename}')
        ndat_cal = self.xffts_cal.getSize()
        calRpos = np.array([x == ('CAL', 'R') for x in self.xffts_cal.scanPattern])
        calSpos = np.array([x == ('CAL', 'S') for x in self.xffts_cal.scanPattern])
        calRspec = []
        calRint  = 0
        calSspec = []
        calSint  =0
        print('   collecting R-pos data ...')
        for ical in tqdm(np.arange(ndat_cal)[calRpos]) :
            calTint = self.xffts_cal.binHead(ical)['integTime']
            calSpec = binnedArray2(self.xffts_cal.readData(ical), binning)
            calRspec.append(calSpec)
            calRint += calTint
        print('   collecting S-pos data ...')
        for ical in tqdm(np.arange(ndat_cal)[calSpos]) :
            calTint = self.xffts_cal.binHead(ical)['integTime']
            calSpec = binnedArray2(self.xffts_cal.readData(ical), binning)
            calSspec.append(calSpec)
            calSint += calTint

        rspec = np.sum( calRspec, axis=0)/calRint
        sspec = np.sum( calSspec, axis=0)/calSint

        Tsys = (273.+self.Temperature_cal) * sspec / (rspec - sspec)

        print('   measured Y-factor (all spw average) is {yf:1.3}'.format(yf = np.log10(np.median(rspec/sspec))*10))
        print('   measured Tsys     (all spw average) is {tsys:3.2f}K'.format(tsys = np.median(Tsys)))

        self.rspec = rspec
        self.sspec = sspec
        self.Tsys = Tsys
        self.Tsys_time = np.median(self.betime_cal.astype('float'))

    def getOnOffSpecs(self, binning=1, interp='linear',noRefCal=True):
        tant_on = self.t_ant.astype(float)[self.antonpos]
        tbe_on  = self.betime.astype(float)[self.beonpos]
        azon_raw = self.az_map[self.antonpos] * 3600 * 180/pi     #rad -> sec
        elon_raw = self.el_map[self.antonpos] * 3600 * 180/pi

        azon_intrp  = interp1d(tant_on - tant_on[0], azon_raw, fill_value='extrapolate')(tbe_on - tant_on[0])
        elon_intrp  = interp1d(tant_on - tant_on[0], elon_raw, fill_value='extrapolate')(tbe_on - tant_on[0])

        azon_raw_sky = self.az_sky[self.antonpos] * 3600 * 180/pi     #rad -> sec
        elon_raw_sky = self.el_sky[self.antonpos] * 3600 * 180/pi

        azon_intrp_sky = interp1d(tant_on - tant_on[0], azon_raw_sky, fill_value='extrapolate')(tbe_on - tant_on[0])
        elon_intrp_sky = interp1d(tant_on - tant_on[0], elon_raw_sky, fill_value='extrapolate')(tbe_on - tant_on[0])


        #### create off-pos spectra #####


        # off pos

        print('### integrating off position spectra ... ')
        offscans =  np.unique(np.array(self.xffts.scanIds)[self.beoffpos] )
        offtime = [None] * len(offscans)
        offspecs= [None] * len(offscans)

        if len(offscans)>0:
            for ix in tqdm(range(len(offscans))) :
                scanid = offscans[ix]
                #time
                avgtime = np.average(self.betime[ (np.array(self.xffts.scanIds)==scanid) * self.ntpsync_work].astype(float))
                offtime[ix] = avgtime

                #spec
                timesum = 0
                specsum = 0
                for idat in np.arange(self.nbedat)[(np.array(self.xffts.scanIds)==scanid) * self.ntpsync_work] :
                    theData = self.xffts.readData(idat)
                    tint    = self.integtime[idat].astype(float)
                    timesum += tint
                    specsum += binnedArray2(theData,binning)  # * tint

                offspecs[ix] = specsum/timesum
        elif noRefCal:
            print(' No off data => Sky level will be estimated from time-stream data.')
        else:
            print(' No off data => skiped')

        # on pos

        onspecIx = np.arange(len(self.beonpos))[self.beonpos]
        ontime = self.betime.astype(float)[self.beonpos]

        print('### collecting on-position spectra ... ')
        onspecs = [None] * len(onspecIx)
        for i in tqdm(range(len(onspecIx))) :
#            dat = self.xffts.readSpec0(onspecIx[i])
#            onspecs[i]= binnedArray2(dat['data'],binning)/dat['integTime']
            dat = self.xffts.readData(onspecIx[i])
            onspecs[i]= binnedArray2(dat,binning)/self.integtime[onspecIx[i]].astype(float)

        if len(offscans)>0:
            print('### interpolating off-position spectra ... ')
            if interp=='linear' :
                offspec_interp = interp1d( offtime - offtime[0], offspecs, axis=0, fill_value='extrapolate') (ontime - offtime[0])
            elif interp=='nearest' :
                offspec_interp = np.array([  offspecs[np.abs( offtime - ontime[i]).argmin()] for i in range(len(ontime))])
            else :
                raise Exception(f'unsupported interpolation {interp}')

            offspec_interp[offspec_interp==0.] = np.nan
            onoffspecs_raw = (onspecs - offspec_interp)/offspec_interp * self.Tsys

            skyzeroflag = np.isnan(onoffspecs_raw).astype('int')
            onoffspecs_raw[np.isnan(onoffspecs_raw)] = 0.

            onoffspecs = np.array([ [ispec[j] - astropy.stats.biweight_location(ispec[j]) for j in range(self.xffts.narray)] for ispec in onoffspecs_raw ])

            self.onoffspecs = onoffspecs
            self.ontime = ontime
            self.skyzeroflag = skyzeroflag
            self.onoffspecs_raw = onoffspecs_raw
            self.offspecs = offspecs

        elif noRefCal:
            print(' Sky level is estimated after masking the signal from the source.')
            offspec_interp = np.zeros_like(onspecs)
            #offspecs_filtered = np.ma.zeros(np.array(onspecs).shape)

            for i in tqdm(range(np.array(onspecs).shape[1])):

                model_sky = models.Polynomial1D(1)
                model_sky.c0 = 1.0e7
                model_sky.c1 = 0.0
                pfit = fitting.LinearLSQFitter()
                opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=15, sigma=2.0)
                fitted_sky,filtered_data = opfit(model_sky, ontime-ontime[0], np.mean(np.array(onspecs)[:,i,:],axis=1))

                if interp=='linear' :
                    offspec_mean = interp1d( ontime[~filtered_data.astype('bool')] - ontime[0], np.mean(np.array(onspecs)[:,i,:],axis=1)[~filtered_data.astype('bool')],fill_value='extrapolate') (ontime - ontime[0])
                elif interp=='nearest' :
                    offspec_mean = np.array([  np.mean(np.array(onspecs)[:,i,:],axis=1)[~filtered_data.astype('bool')][np.abs( ontime[~filtered_data.astype('bool')] - ontime[i]).argmin()] for ii in range(len(ontime))])
                elif interp=='fit':
                    offspec_mean = fitted_sky(ontime-ontime[0])
                else :
                    raise Exception(f'unsupported interpolation {interp}')


                offspec_interp[:,i,:] = np.array([np.full(np.array(onspecs).shape[2],offspec_mean[ii]) for ii in range(ontime.shape[0])])

            offspec_interp[offspec_interp==0.] = np.nan
            onoffspecs_raw = (onspecs - offspec_interp)/offspec_interp * self.Tsys
            onoffspecs_cont = (np.mean(np.array(onspecs),axis=2) - np.mean(offspec_interp,axis=2))/np.mean(offspec_interp,axis=2) * np.median(self.Tsys,axis=1)

            skyzeroflag = np.isnan(onoffspecs_raw).astype('int')
            onoffspecs_raw[np.isnan(onoffspecs_raw)] = 0.

            onoffspecs = np.array([ [ispec[j] - astropy.stats.biweight_location(ispec[j]) for j in range(self.xffts.narray)] for ispec in onoffspecs_raw ])

            #self.offspecs_filtered = offspecs_filtered
            self.onoffspecs_cont = onoffspecs_cont
            self.onoffspecs = onoffspecs
            self.offspec_interp = offspec_interp
            self.skyzeroflag = skyzeroflag
            self.onoffspecs_raw = onoffspecs_raw
            self.offspecs = offspecs


        else:
            print(' No off data => skiped')

        self.onspecs = onspecs
        self.ontime = ontime
        self.lsbfreq       = binnedArray(self.flsb, binning)
        self.usbfreq       = binnedArray(self.fusb, binning)
        self.azon_intrp = azon_intrp
        self.elon_intrp = elon_intrp
        self.azon_intrp_sky = azon_intrp_sky
        self.elon_intrp_sky = elon_intrp_sky
        self.azon_raw   = azon_raw
        self.elon_raw   = elon_raw

    def mkMap(self, useSBForMap='lsb', usePolForMap=['pol1'], size=[None, None], grid=[3.,3.], smooth=[3.,3.], chrangeMap=[None, None], searchRadius=2. , mode='int',cont=False):

        if cont:
            onoffspecs = np.zeros([self.onoffspecs_cont.shape[0],self.onoffspecs_cont.shape[1],1])
            onoffspecs[:,:,0] = self.onoffspecs_cont
            chrangeMap = [0,1]
        else:
            onoffspecs = self.onoffspecs
        azon_raw = self.az_map[self.antonpos] * 3600 * 180/np.pi     #rad -> sec
        elon_raw = self.el_map[self.antonpos] * 3600 * 180/np.pi
        tant_on = self.t_ant.astype(float)[self.antonpos]
        self.smoothFWHM = np.array(smooth)
        smoothSTD = self.smoothFWHM/2.35

        mpx = self.azon_intrp
        mpy = self.elon_intrp

        if size[0] == None :
            xmax, xmin = np.floor(mpx.max()/grid[0]) * grid[0],  np.floor(mpx.min()/grid[0])*grid[0]
        else :
            xmax, xmin = size[0]/2., -size[0]/2.

        if size[1] == None :
            ymax, ymin = np.floor(mpy.max()/grid[1]) * grid[1],  np.floor(mpy.min()/grid[1])*grid[1]
        else :
            ymax, ymin = size[1]/2., -size[1]/2.

        usespwForMap = [ i for i in range(self.nIF) if (self.sbdef[i] == useSBForMap ) and ( self.poldef[i] in usePolForMap) ]

        if chrangeMap[0] == None:
            chrangeMap[0] = 0
        if chrangeMap[1] == None:
            if useSBForMap=='lsb':
                chrangeMap[1] = len(self.lsbfreq)
            elif useSBForMap=='usb':
                chrangeMap[1] = len(self.usbfreq)

        print(chrangeMap)

        nx = int(abs(xmax-xmin)/grid[0]) + 1
        ny = int(abs(ymax-ymin)/grid[1]) + 1
        nz = chrangeMap[1] - chrangeMap[0]


        ###########
        # create axis

        xaxis = np.linspace( xmin, xmax, num=nx, endpoint=True)
        yaxis = np.linspace( ymin, ymax, num=ny, endpoint=True)

        ###########
        #

        imagemap  = np.zeros((nx,ny,nz))
        weightmap = np.zeros((nx,ny))

        ###########
        # gaussian
        #print(self.Ta.copy().values.shape)
        gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
        obs_sp = np.mean(onoffspecs[:,usespwForMap, chrangeMap[0]:chrangeMap[1]],axis=1)
        kpix    = float(self.smoothFWHM[0])/grid[0], float(self.smoothFWHM[1])/grid[1]
        search  = int(kpix[0]*searchRadius), int(kpix[1]*searchRadius)
        print(kpix, search)

        ##########

        for i in tqdm(range(obs_sp.shape[0])):
            xix = np.abs( mpx[i] - xaxis).argmin()
            yix = np.abs( mpy[i] - yaxis).argmin()

            xr = max([xix-search[0],0]), min([xix+search[0],nx-1])
            yr = max([yix-search[1],0]), min([yix+search[1],ny-1])

            wtx = gsbeam( xaxis[xr[0]:xr[1]] - mpx[i], smoothSTD[0])[:,np.newaxis]
            wty = gsbeam( yaxis[yr[0]:yr[1]] - mpy[i], smoothSTD[1])[:,np.newaxis]

            wtmap_i = np.dot(wtx, wty.T)
            weightmap[xr[0]:xr[1],yr[0]:yr[1]] += wtmap_i

            imagemap_i = wtmap_i[:,:,np.newaxis] * obs_sp[i]
            imagemap[xr[0]:xr[1],yr[0]:yr[1]] += imagemap_i

        masked = (weightmap==0.)

        weightmap[masked] = 0.1

        self.imCube  = imagemap/weightmap[:,:,np.newaxis]
        if mode == 'int' :   #integrated intensity
            self.imObj   = np.sum(self.imCube,axis=2)
        elif mode == 'peak' :  #peak intensity
            self.imObj   = np.max    (self.imCube,axis=2)
        else :
            raise Exception(f'unsupported mode {mode}')

#        self.imObj[masked] = None
        self.obspk = np.max(self.imObj)
        obspix= np.where( self.imObj == self.obspk )
        self.obspos= xaxis[obspix[0][0]], yaxis[obspix[1][0]]
        self.obsrms= astropy.stats.biweight_scale( self.imObj[masked==False] )
        self.obssky= astropy.stats.biweight_location( self.imObj[masked==False] )
        self.xaxis=xaxis
        self.yaxis=yaxis
        self.mask = masked
        self.weightmap = weightmap
        self.mkMapMode = mode


        print("   pos        = {f:}".format(f=self.obspos))
        print("   pk         = {f:.2f}".format(f=self.obspk))
        print("   source/rms = {f:.2f}".format(f=self.obspk/self.obsrms))
        print("   source/sky = {f:.2f}".format(f=self.obspk/self.obssky))

        #fit = plt.figure(figsize=(12,5))
        '''
        plt.subplot(3,2,3)
        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.title('{name}'.format(name={'int':'Integrated Intensity', 'peak':'peak intensity'}[mode]))
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.subplot(3,2,4)
        plt.imshow(np.rot90(weightmap)[:,::-1], extent=(xmax, xmin, ymin, ymax), cmap='jet')
        plt.title('weight')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.savefig('./test.png')
        '''

    def mkFit(self) :
        xmax,xmin = np.max(self.xaxis), np.min(self.xaxis)
        ymax,ymin = np.max(self.yaxis), np.min(self.yaxis)

        p0=(self.obspos[0], self.obspos[1],
            10/2.35, 10/2.35,
            self.obspk , self.obssky )

        popt, trash = scipy.optimize.curve_fit( modelBeam, np.meshgrid(self.xaxis,self.yaxis), self.imObj.flatten() , p0=p0, maxfev=100000)
        modelim = modelBeam(np.meshgrid(self.xaxis, self.yaxis), *popt).reshape(self.imObj.shape)
        self.fitloc  = np.array([popt[0], popt[1]])
        self.fitbeam = np.array([popt[2], popt[3]])*2.35
        self.fitbeam_corr = np.sqrt( self.fitbeam**2 - self.smoothFWHM**2)
        self.fitpeak = (popt[4] - popt[5])
        self.fitpeak_err = np.sqrt(np.sqrt(np.diag(trash))[4]**2 + np.sqrt(np.diag(trash))[5]**2)
        self.modelim = modelim

        print('### beam fitting result :')
        print('   location       : {x:1.3f}",  {y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1]))
        print('   FWHM(raw)      : {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam[0], y=self.fitbeam[1]))
        print('   FWHM(corrected): {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1]))

    def searchLine(self, nTimeAvg=8, vrange =[ -200, 200], sig_thres=3.,useFreq=129.363,useSBForSpec = 'lsb',usePolForSpec=['pol1']) :

        ##
        print('### Searching line(s) ... ')
        if useSBForSpec == 'lsb':
            freq = self.lsbfreq
        elif useSBForSpec == 'usb':
            freq = self.usbfreq

        if useFreq < freq.min() or useFreq > freq.max():
            useFreq=freq.mean()
            print(' [WARN] useFreq='+str(useFreq)+'GHz is out of the frequency range.')
            print(' [WARN] useFreq='+str(freq.maen())+'GHz (band center) will be used.')

        vel = ( 1 - freq/useFreq) * 2.99792458e5

        if vrange ==[ None, None] :
            fitch = [ 0, len(freq)-1 ]

        else :
            fitch = sorted([ np.abs( vel - vrange[0]).argmin(), np.abs( vel - vrange[1]).argmin() ])

        usespwForSpec = [ i for i in range(self.nIF) if (self.sbdef[i] == useSBForSpec ) and ( self.poldef[i] in usePolForSpec) ]
        obs_sp = np.mean(self.onoffspecs[:,usespwForSpec,],axis=1)

        specTimeAvg = binnedArray(obs_sp, nTimeAvg)
        specMax = specTimeAvg[np.argmax(np.max(specTimeAvg[:,fitch[0]:fitch[1]],axis=1))]

        filtered_data = sigma_clip(specMax[fitch[0]:fitch[1]], cenfunc='mean',sigma=sig_thres, masked=True,maxiters=15).mask
        bl = np.mean(specMax[fitch[0]:fitch[1]][~filtered_data.astype('bool')])
        chRange_vel = checkspec(filtered_data,startch=np.argmax(specMax[fitch[0]:fitch[1]]))
        chRange = [chRange_vel[0]+fitch[0],chRange_vel[1]+fitch[0]]

        self.vel = vel
        self.cenfreq = useFreq
        self.spec_vel = specMax
        self.chRange = chRange
        self.bl = bl

    def efficiency(self,D = 50.,useSBForEff='lsb',usePolForEff='pol1'):
        if useSBForEff == 'lsb':
            self.useFreq_cont = np.mean(self.lsbfreq)
        elif useSBForEff == 'usb':
            self.useFreq_cont = np.mean(self.usbfreq)

        HPBW_x = self.fitbeam_corr[0]
        HPBW_y = self.fitbeam_corr[1]
        lobs = 3.0e8/self.useFreq_cont/1000.
        #Tpeak = 3.84641
        #amp = 4.13259

        # Brightness temperature of Uranus (Griffin & Orton 1993)
        a0 = -795.694
        a1 = 845.179
        a2 = -288.946
        a3 = 35.200
        Turanus = a0 + a1*np.log10(lobs) + a2*(np.log10(lobs))**2 + a3*(np.log10(lobs))**3 # [K]
        Turanus = Turanus*0.931 # Sayers+12

        # main beam efficiency
        theta_eq = 3.72 # [arcsec]
        theta_pol = 3.64 # [arcsec]

        #bff = 1 - np.exp(-np.log(2) * theta_eq * theta_pol / fwhm_x / fwhm_y)
        #eta = Tpeak/Turanus/bff
        bff = 1 - np.exp(-np.log(2) * theta_eq * theta_pol / HPBW_x/HPBW_y)
        eta = self.fitpeak/Turanus/bff
        eta_err = self.fitpeak_err/Turanus/bff

        # aperture efficiency
        etaa = eta/1.2 # assuming -12 dB edge taper
        etaa_err = eta_err/1.2

        # Jy/K
        JpK = HPBW_x*HPBW_y/13.6/(lobs*1e-3)**2

        self.eta = eta
        self.etaa = etaa
        self.JpK = JpK
        self.eta_err = eta_err
        self.etaa_err = etaa_err
        self.Uranus_K = Turanus
        self.Uranus_Jy = Turanus*JpK*bff
        self.bff = bff

    def createPsData(self,Dir='PsData/'):
        os.system('mkdir -p '+globLogDir)
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum))
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum)+'/'+Dir)
        npyfile=globLogDir+'/'+str(self.obsnum)+'/'+Dir+'PsData.'+str(self.obsnum)+'.chbin'+str(self.binning)+'.npy'
        np.save(npyfile,self.onoffspecs_raw)

        f = open(npyfile+'.readme.txt','w')
        f.write('0: time [id]'+'\n')
        f.write('1: spw [0-3] # 0=Bpol LSB, 1=Bpol USB, 2=Apol LSB, 3=Apol USB'+'\n')
        f.write('2: freq [GHz]'+'\n')
        f.close()

        np.save(npyfile.replace('.npy','.lsbfreq.npy'),self.lsbfreq.data)
        np.save(npyfile.replace('.npy','.usbfreq.npy'),self.usbfreq.data)


    def Pipeline(self,binning=8,interp='linear',noRefCal=False):
        self.loadlmttpm()
        self.loadxffts()
        self.TsysCal(binning=binning)
        self.getOnOffSpecs(binning=binning,interp=interp,noRefCal=noRefCal)
        self.binning = binning

    def PsQlook(self,AvgScan=10,rmsflag_thrs=1.2,highz=False):

        os.system('mkdir -p '+globLogDir)
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum))
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum)+'/PsQlook')

        onoffspecs_raw_avg = binnedArray(self.onoffspecs_raw,AvgScan)
        onoffspecs_bl_avg = np.zeros_like(onoffspecs_raw_avg)

        Tsys_value = np.median(self.Tsys,axis=1)

        for nn in range(self.nIF):
            print('### Ps Qlook for '+self.sbdef[nn]+' '+self.poldef[nn]+'...')
            numflag=0
            rmsflag = ~np.zeros(onoffspecs_raw_avg.shape[0],dtype='bool')

            for i in tqdm(range(onoffspecs_raw_avg.shape[0])):
                if self.sbdef[nn]=='lsb':
                    freq = self.lsbfreq
                elif self.sbdef[nn]=='usb':
                    freq = self.usbfreq

                model_bl = models.Polynomial1D(1)
                model_bl.c0 = 0.0
                model_bl.c1 = 0.0001
                pfit = fitting.LinearLSQFitter()
                opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=15, sigma=3.0)
                fitted_bl,filtered_data = opfit(model_bl, freq, onoffspecs_raw_avg[i][nn])
                onoffspecs_bl_avg[i][nn] = onoffspecs_raw_avg[i][nn]-fitted_bl(freq)
                rms = np.std( onoffspecs_bl_avg[i][nn][~filtered_data.astype('bool')] )

                if self.integtime.astype('float').min()/1.0e6 < 0.5:
                    rms_expect = RMSexpect(Tsys_value[nn],np.abs(np.mean(np.diff(freq))),tint=AvgScan/5.)
                else:
                    rms_expect = RMSexpect(Tsys_value[nn],np.abs(np.mean(np.diff(freq))),tint=AvgScan)

                if rms > rmsflag_thrs*rms_expect:
                    rmsflag[i] = False
                    numflag=numflag+1

            print('   '+str(numflag)+'/'+str(onoffspecs_raw_avg.shape[0])+' AvgScans flagged.')

            spec = np.mean(onoffspecs_bl_avg[:,nn,:][rmsflag],axis=0)

            plt.close()
            plt.title(self.srcname+' ('+self.sbdef[nn]+','+['Bpol','Bpol','Apol','Apol'][nn]+') chbin: '+str(self.binning)+' Tsys: {x:1.2f}K'.format(x=Tsys_value[nn]))
            if highz:
                plt.fill_between(freq,np.zeros_like(freq),spec,step='mid',facecolor='b')
            else:
                plt.step(freq,spec,where='mid')
            plt.xlabel('freq [GHz]')
            plt.ylabel('Ta* [K]')
            plt.savefig(globLogDir+'/'+str(self.obsnum)+'/PsQlook/PsQlook.'+str(self.obsnum)+'.chbin'+str(self.binning)+'.'+str(nn+1).zfill(2)+'.png')
            np.save(globLogDir+'/'+str(self.obsnum)+'/PsQlook/PsQlook.'+str(self.obsnum)+'.chbin'+str(self.binning)+'.'+str(nn+1).zfill(2)+'.npy',np.vstack([freq,spec]).T)
            plt.close()

    def LinePointingQlook(self,showVelRange=[-150,150]):
        self.searchLine()
        self.mkMap(chrangeMap=[self.chRange[0],self.chRange[1]])
        self.mkFit()

        plt.close()
        plt.rcParams["font.size"] = 7
        plt.figure(figsize=[6,6.5])
        plt.subplots_adjust(wspace=0.3, hspace=0.5,top=0.9,bottom=0.1)
        plt.subplot(3,1,1)

        titleinfo = '['+self.srcname+'] offset: {x:1.3f}",{y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1])+' FWHM(corr): {x:1.3f}"x{y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1])
        plt.title(titleinfo)

        plt.step( self.vel, self.spec_vel-self.bl, where='mid')
        plt.xlim(showVelRange[0],showVelRange[1])
        yrange = (min(self.spec_vel-self.bl), max(self.spec_vel-self.bl))
        plt.ylim( (yrange[0]*1.1 - yrange[1]*0.1, yrange[1]*1.1 - yrange[0]*0.1 ))
        plt.xlabel('V_topo (km/s)')
        plt.grid()

        plt.fill_between( self.vel[self.chRange[0]:self.chRange[1]], self.spec_vel[self.chRange[0]:self.chRange[1]]-self.bl, color='r', step='mid')

        xmax,xmin = np.max(self.xaxis), np.min(self.xaxis)
        ymax,ymin = np.max(self.yaxis), np.min(self.yaxis)

        plt.subplot(3,2,3)
        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.title('{name}'.format(name={'int':'Integrated Intensity', 'peak':'peak intensity'}[self.mkMapMode]))
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.subplot(3,2,4)
        plt.imshow(np.rot90(self.weightmap)[:,::-1], extent=(xmax, xmin, ymin, ymax), cmap='jet')
        plt.title('weight')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        plt.subplot(3,2,5)
        plt.title('model')
        plt.imshow(np.rot90(self.modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.subplot(3,2,6)
        plt.title('residual')
        plt.imshow(np.rot90(self.imObj - self.modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        os.system('mkdir -p '+globLogDir)
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum))
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum)+'/LinePointing')
        plt.savefig(globLogDir+'/'+str(self.obsnum)+'/LinePointing/LinePointing.'+str(self.obsnum)+'.png')
        plt.close()

        SourceAZEL_index = np.argmin(np.sqrt(self.az_map**2+self.el_map**2))
        AZ = np.rad2deg(self.source_az)[SourceAZEL_index]
        EL = np.rad2deg(self.source_el)[SourceAZEL_index]
        dAZuser = np.rad2deg(self.dAZuser)*3600.
        dELuser = np.rad2deg(self.dELuser)*3600.
        dAZ_model = np.rad2deg(self.dAZModel)*3600.
        dEL_model = np.rad2deg(self.dELModel)*3600.
        dAZ_M2 = np.rad2deg(self.dAZM2Cor )*3600.
        dEL_M2 = np.rad2deg(self.dELM2Cor)*3600.

        f = open(globLogDir+'/'+str(self.obsnum)+'/LinePointing/LinePointing.'+str(self.obsnum)+'.log','w')
        f.write('#(0)obsnum (1)AZ (2)EL (3)dAZ(user) (3)dAZ(user) (4)dEL(user) (5)dAZ(model) (6)dEL(model) (7)dAZ(M2) (8)dEL(M2) (9)Temp[degC] (10)Tau@220GHz (11)WindSpeed1 (12)WindDir1 (13)WindSpeed2 (14)WindDir2'+'\n')
        f.write('# Source Name: '+self.srcname+'\n')
        f.write('# AZEL unit: arcsec'+'\n')
        f.write(str(AZ)+' ')
        f.write(str(EL)+' ')
        f.write(str(dAZ_model)+' ')
        f.write(str(dEL_model)+' ')
        f.write(str(dAZ_M2)+' ')
        f.write(str(dEL_M2)+' ')
        f.write(str(self.Temperature)+' ')
        f.write(str(self.tau)+' ')
        f.write(str(self.WindSpeed1)+' ')
        f.write(str(self.WindDir1)+' ')
        f.write(str(self.WindSpeed2)+' ')
        f.write(str(self.WindDir2)+'\n')
        f.close()

        self.createPsData(Dir='LinePointing/')

    def ContPointingQlook(self):
        #for i in range(self.nIF):
        for i in range(self.nIF):
            print('')
            print('### '+self.sbdef[i]+' '+self.poldef[i]+'... ')
            self.mkMap(useSBForMap=self.sbdef[i], usePolForMap=[self.poldef[i]],cont=True)
            self.mkFit()

            plt.close()
            plt.rcParams["font.size"] = 7
            plt.figure(figsize=[6,8.7])
            plt.subplots_adjust(wspace=0.3, hspace=0.5,top=0.9,bottom=0.1)
            plt.subplot(4,1,1)

            titleinfo = '['+self.srcname+']'+' ('+self.sbdef[i]+','+['Bpol','Bpol','Apol','Apol'][i]+')'+' offset: {x:1.3f}",{y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1])+' FWHM(corr): {x:1.3f}"x{y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1])
            plt.title(titleinfo)

            plt.plot((self.ontime-self.ontime[0])/1.0e6,np.mean(np.array(self.onspecs)[:,i,],axis=1)/1.0e7,'ro-',markersize=2.,label='Signal')
            plt.plot((self.ontime-self.ontime[0])/1.0e6,np.mean(self.offspec_interp[:,i,:],axis=1)/1.0e7,'bo-',markersize=1.5,label='Estimated sky')

            plt.xlabel('time (sec)')
            plt.ylabel('XFFTS output * 1.0e7')
            plt.grid()
            plt.legend()

            plt.subplot(4,1,2)

            titleinfo = 'Peak(Map): {x:1.3f}'.format(x=self.fitpeak)+'K Peak(ts): {x:1.3f}'.format(x=np.max(self.onoffspecs_cont[:,i]))+'K'
            if self.srcname=='Uranus':
                self.efficiency(useSBForEff=self.sbdef[i], usePolForEff=[self.poldef[i]])
                titleinfo = titleinfo + '  eta_a(Map)={x:1.3f}'.format(x=self.etaa)+' eta_a(ts)={x:1.3f}'.format(x=self.etaa*np.max(self.onoffspecs_cont[:,i])/self.fitpeak)+' at {x:3.1f}GHz'.format(x=self.useFreq_cont)
            plt.title(titleinfo)

            plt.plot((self.ontime-self.ontime[0])/1.0e6,self.onoffspecs_cont[:,i],'ro-',markersize=2.,label='Signal')

            plt.xlabel('time (sec)')
            plt.ylabel('Ta* [K]')
            plt.grid()

            xmax,xmin = np.max(self.xaxis), np.min(self.xaxis)
            ymax,ymin = np.max(self.yaxis), np.min(self.yaxis)

            plt.subplot(4,2,5)
            plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
            plt.title('{name}'.format(name={'int':'Integrated Intensity', 'peak':'peak intensity'}[self.mkMapMode]))
            plt.xlabel('dAz (arcsec)')
            plt.ylabel('dEl (arcsec)')
            plt.subplot(4,2,6)
            plt.imshow(np.rot90(self.weightmap)[:,::-1], extent=(xmax, xmin, ymin, ymax), cmap='jet')
            plt.title('weight')
            plt.xlabel('dAz (arcsec)')
            plt.ylabel('dEl (arcsec)')

            plt.subplot(4,2,7)
            plt.title('model')
            plt.imshow(np.rot90(self.modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
            plt.xlabel('dAz (arcsec)')
            plt.ylabel('dEl (arcsec)')
            plt.subplot(4,2,8)
            plt.title('residual')
            plt.imshow(np.rot90(self.imObj - self.modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
            plt.xlabel('dAz (arcsec)')
            plt.ylabel('dEl (arcsec)')

            os.system('mkdir -p '+globLogDir)
            os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum))
            os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum)+'/ContPointing')
            plt.savefig(globLogDir+'/'+str(self.obsnum)+'/ContPointing/ContPointing.'+str(self.obsnum)+'.'+str(i+1).zfill(2)+'.png')
            plt.close()

            SourceAZEL_index = np.argmin(np.sqrt(self.az_map**2+self.el_map**2))
            AZ = np.rad2deg(self.source_az)[SourceAZEL_index]
            EL = np.rad2deg(self.source_el)[SourceAZEL_index]
            dAZuser = np.rad2deg(self.dAZuser)*3600.
            dELuser = np.rad2deg(self.dELuser)*3600.
            dAZ_model = np.rad2deg(self.dAZModel)*3600.
            dEL_model = np.rad2deg(self.dELModel)*3600.
            dAZ_M2 = np.rad2deg(self.dAZM2Cor )*3600.
            dEL_M2 = np.rad2deg(self.dELM2Cor)*3600.

            f = open(globLogDir+'/'+str(self.obsnum)+'/ContPointing/ContPointing.'+str(self.obsnum)+'.'+str(i+1).zfill(2)+'.log','w')
            f.write('#(0)obsnum (1)AZ (2)EL (3)dAZ(user) (3)dAZ(user) (4)dEL(user) (5)dAZ(model) (6)dEL(model) (7)dAZ(M2) (8)dEL(M2) (9)Temp[degC] (10)Tau@220GHz (11)WindSpeed1 (12)WindDir1 (13)WindSpeed2 (14)WindDir2'+'\n')
            f.write('# Source Name: '+self.srcname+'\n')
            f.write('# AZEL unit: arcsec'+'\n')
            f.write(str(AZ)+' ')
            f.write(str(EL)+' ')
            f.write(str(dAZ_model)+' ')
            f.write(str(dEL_model)+' ')
            f.write(str(dAZ_M2)+' ')
            f.write(str(dEL_M2)+' ')
            f.write(str(self.Temperature)+' ')
            f.write(str(self.tau)+' ')
            f.write(str(self.WindSpeed1)+' ')
            f.write(str(self.WindDir1)+' ')
            f.write(str(self.WindSpeed2)+' ')
            f.write(str(self.WindDir2)+'\n')
            f.close()

            f = open(globLogDir+'/'+str(self.obsnum)+'/ContPointing/ContPointing.'+str(self.obsnum)+'.eff.'+str(i+1).zfill(2)+'.log','w')
            f.write('#(0)obsnum (1)Freq [GHz] (2)Peak(Map)[K] (3)eta_a(Map) (4)Peak(Map)[K] (5)eta_a(Map)'+'\n')
            f.write(str(self.obsnum)+' ')
            f.write(str(self.useFreq_cont)+' ')
            f.write(str(self.fitpeak)+' ')
            f.write(str(self.etaa)+' ')
            f.write(str(np.max(self.onoffspecs_cont[:,i]))+' ')
            f.write(str(self.etaa*np.max(self.onoffspecs_cont[:,i])/self.fitpeak)+'\n')
            f.close()

            self.createPsData(Dir='ContPointing/')

    def PointingConv(self,dAZ=0.,dEL=0.):
        print('### Convert AZEL to RADEC...')
        # dAZ, dEL in the arcsec unit

        direction = np.zeros((self.ontime.shape[0],2),dtype='float64')
        time = np.zeros(self.ontime.shape[0],dtype='float64')

        for i in tqdm(range(self.ontime.shape[0])):
            onTime = Time(self.ontime[i]/1.0e6,format='unix',scale='utc')
            ra,dec = AZEL2RADEC((self.azon_intrp_sky[i]-dAZ)/3600.,(self.elon_intrp_sky[i]-dEL)/3600.,onTime)
            direction[i][0] = ra
            direction[i][1] = dec
            time[i] = onTime.mjd*24.*60.*60.

        self.direction = direction
        self.time = time
        self.dAZ = dAZ
        self.dEL = dEL

    def createMS2(self,dAZ=0.,dEL=0.):
        self.PointingConv(dAZ=dAZ,dEL=dEL)
        specdata = np.zeros([self.time.shape[0],2,2,self.lsbfreq.shape[0]])
        specdata[:,0,0,:] = self.onoffspecs_raw[:,0,:]
        specdata[:,0,1,:] = self.onoffspecs_raw[:,1,:]
        specdata[:,1,0,:] = self.onoffspecs_raw[:,2,:]
        specdata[:,1,1,:] = self.onoffspecs_raw[:,3,:]

        Tsys = np.zeros([2,2,self.lsbfreq.shape[0]])
        Tsys[0,0,:] = self.Tsys[0,:]
        Tsys[0,1,:] = self.Tsys[1,:]
        Tsys[1,0,:] = self.Tsys[2,:]
        Tsys[1,1,:] = self.Tsys[3,:]

        Tsys_time = Time(self.Tsys_time/1.0e6,format='unix',scale='utc').mjd*24.*60.*60.
        state_id = np.full(self.time.shape[0],1).astype('int')
        freq = np.vstack([self.lsbfreq,self.usbfreq])*1.0e9

        os.system('mkdir -p '+globLogDir)
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum))
        os.system('mkdir -p '+globLogDir+'/'+str(self.obsnum)+'/MS2')

        if self.dAZ==0. and self.dEL==0.:
            MS2name = globLogDir+'/'+str(self.obsnum)+'/MS2/MS2.'+str(self.obsnum)+'.chbin'+str(self.binning)+'.ms'
        else:
            MS2name = globLogDir+'/'+str(self.obsnum)+'/MS2/MS2.'+str(self.obsnum)+'.chbin'+str(self.binning)+'.dAZ_'+str(self.dAZ)+'.EL_'+str(self.dEL)+'.ms'

        ms2.createMS2(MS2name,specdata,self.time,np.rad2deg(self.source_ra),np.rad2deg(self.source_dec),self.sysvel,self.srcname,self.Pid,self.Observer,self.direction,freq,Tsys,Tsys_time,state_id)

###
