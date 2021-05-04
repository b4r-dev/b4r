import glob
import re
import netCDF4 as nc
import b4rpipe.specfile4 as sp
import numpy as np
import math as ma


def getEnvironment():
    try:
        env = get_ipython().__class__.__name__
        if env == 'ZMQInteractiveShell':
            return 'Jupyter'
        elif env == 'TerminalInteractiveShell':
            return 'IPython'
        else:
            return 'OtherShell'
    except NameError:
        return 'Interpreter'



isnoteBook = getEnvironment() == 'Jupyter'
if isnoteBook:
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm


def radec2str(ra,dec, inunit='rad') :
#    print(inunit, ra,dec)

    if inunit == 'rad' :
        radeg  = ra  * 180/ma.pi
        decdeg = dec * 180/ma.pi
    else :
        radeg = ra
        decdeg = dec

    rasig = 1 if radeg > 0 else -1
    rahabs = radeg / 15. * rasig
    rah   = int(rahabs)
    ram   = int((rahabs - rah)*60)
    ras   = (rahabs - rah - ram/60.)*3600

    decsig = 1 if decdeg > 0 else -1
    decabs = decdeg  * decsig
    decd   = int(decabs)
    decm   = int((decabs - decd)*60)
    decs   = (decabs - decd - decm/60.)*3600

    return ( '{sig}{h:02d}h{m:02d}m{s:02.2f}s'.format(sig=' ' if rasig == 1 else '-',  h=rah, m=ram, s=ras),
             '{sig}{h:02d}d{m:02d}m{s:02.2f}s'.format(sig=' ' if rasig == 1 else '-',  h=decd, m=decm, s=decs) )

def getFileList(basedir, ffilter=None):

#    if not 'basedir' in globals() :
#        basedir = '../'


    lmtdir   = '{d}/lmttpm/'.format(d=basedir)
    xfftsdir = '{d}/xffts/'.format(d=basedir)

    print(f'making file list in = {lmtdir}')

    files = glob.glob('{d}/lmttpm*.nc'.format(d=lmtdir))
    if ffilter != None :
        flt = re.compile(ffilter)
        files_matched = [ x for x in files if flt.search(x) != None ]
        files = files_matched

    files_out = []
    for itr in tqdm(files) :
        try :
            ncf = nc.Dataset(itr)
            if 'Header.B4r.LineFreq' in ncf.variables :
                files_out.append(itr)
            ncf.close()
        except :
            continue

    return files_out

def getSourceList(ffilter=None) :
    files = getFileList(ffilter)

    ############################
    # sources
    ############################

    sources = {}

    for x in files :
        f1 = nc.Dataset(x)
        name = ''.join([  x.decode() for x in f1.variables['Header.Source.SourceName'][:].tolist() if x != b' ' ])
        ra = f1.variables['Header.Source.Ra'][0].tolist()
        dec= f1.variables['Header.Source.Dec'][0].tolist()
        sources[name] =  ( ra,dec)
        f1.close()

    print(sources)

def getFreqList(ffilter=None) :

    files = getFileList(ffilter)

    ############################
    # freqs
    ############################

    obsfreqs = set([])

    for x in files :
        f1 = nc.Dataset(x)
        freq= f1.variables['Header.B4r.LineFreq'][0].tolist()
        obsfreqs.add(freq)
        f1.close()
    print(obsfreqs)


def getNCprops(ncfile) :

    try :
        ncf = nc.Dataset(ncfile)
        srcname = ''.join([  x.decode() for x in ncf.variables['Header.Source.SourceName']         [:].tolist() if x != b' ' ])
        projid = ''.join([  x.decode() for x in ncf.variables['Header.Dcs.ProjectId']         [:].tolist() if x != b' ' ])
        observer = ''.join([  x.decode() for x in ncf.variables['Header.Telescope.Operator']         [:].tolist() if x != b' ' ])
        fmlofile= ''.join([  x.decode() for x in ncf.variables['Header.B4r.FreqModulationFilename'][:].tolist() if x != b' '])
        srcra   = ncf.variables['Header.Source.Ra' ][0].tolist()
        srcdec  = ncf.variables['Header.Source.Dec'][0].tolist()
        radecstr = radec2str(srcra, srcdec);
        tm_st   = ncf.variables['Data.TelescopeBackend.TelTime'][0].tolist()
        tm_ed   = ncf.variables['Data.TelescopeBackend.TelTime'][-1].tolist()
        obsnum  = ncf.variables['Header.Dcs.ObsNum'][0].tolist()
        obsgoal = ''.join([  x.decode() for x in ncf.variables['Header.Dcs.ObsGoal']         [:].tolist() if x != b' ' ])
        freq     = ncf.variables['Header.B4r.LineFreq'][0].tolist()
 #calmode  = ncf.variables['Header.Dcs.CalMode'][0].tolist()
        bufpos   = sorted(np.unique(ncf.variables['Data.TelescopeBackend.BufPos'][:].tolist()))
        obsmode  = 'CALIB' if bufpos == [2,3] else 'OBS' if 0 in bufpos else 'OTHER'
#        azmap  = ncf.variables['Data.TelescopeBackend.TelAzMap'][:]
#        elmap  = ncf.variables['Data.TelescopeBackend.TelElMap'][:]

        ncf.close()

        return {'srcname':srcname, 'srcPos':radecstr, 'tm':[tm_st, tm_ed], 'obsnum':obsnum, 'obsgoal':obsgoal,
                'bufpos':bufpos, 'obsmode':obsmode, 'obsfreq':freq, 'fmlofile':fmlofile, 'observer':observer, 'projid':projid }

    except :
        return None


def getXFFTSbyObsNum(obsnum, basedir) :

    lmtdir   = '{d}/lmttpm/'.format(d=basedir)
    xfftsdir = '{d}/xffts/'.format(d=basedir)

    files = np.unique( [ x[:x.find('.xfftsx')]+'.xfftsx' for x in  glob.glob('{d}/xffts*'.format(d=xfftsdir)) ])

    for x in files :
        try :
            xff = sp.SpecFile(x)
            if xff.binHead(0)['obsNum'] == obsnum :
                return x
        except:
            continue

    return None


def getLMTtpmByObsNum(obsnum, basedir):

    files = getFileList(ffilter=f'_{obsnum:06}_[0-9]{{2}}_[0-9]{{4}}.nc$', basedir=basedir)
    print(files)
    return files[0]
