from ftplib import FTP
import b4rpipe.specfile4 as sp
import sys
import os

class B4Rquery:

    def __init__(self,
                username,
                password,
                hosturl='ngc253.mtk.ioa.s.u-toko.ac.jp',
                b4rDir='/B4R/',
                localBaseDir='./'
                ):

        config = {
            'host':hosturl,
            'user':username,
            'passwd':password,
            }

        self.config = config
        self.ftp = FTP(**config)
        self.localBaseDir = localBaseDir
        self.b4rDir = b4rDir
        self.lmttpm = self.ftp.nlst(b4rDir+'lmttpm/')
        self.xffts = self.ftp.nlst(b4rDir+'xffts/')


    def FindRawdata(self,obsid):
        self.obsid = obsid

        # lmttpm
        result_lmttpm = [s for s in self.mttpm if str(obsid).zfill(6) in s]

        if len(result_lmttpm)==1:
            print('Find one lmttpm file: '+result_lmttpm[0])
            self.lmttpmFile = result_lmttpm[0]
        else:
            print('Error: No or more than one lmttpm files were found', file=sys.stderr)
            sys.exit(1)

        # xffts
        for x in xffts:
            try:
                xff = sp.SpecFile(x)
                if xff.binHead(0)['obsNum'] == obsid:
                    print('Find one xffts file: '+x)
                    self.xfftsFile = x
                    break
            except:
                continue

        if self.xfftsFile==None:
            print('Error: No xffts file was found', file=sys.stderr)
            sys.exit(1)


    def DownloadRawdata(self):

        os.system('mkdir -p '+self.localBaseDir+'lmttpm/')
        os.system('mkdir -p '+self.localBaseDir+'xffts/')

        #lmttpm
        print('Downloading a lmttpm file of obsid: '+str(self.obsid).zfill(6))
        try:
            with FTP(**self.config) as ftp:
                with open(self.lmttpm, 'wb') as fp:
                    ftp.retrbinary('RETR '+elf.localBaseDir+'lmttpm/'+os.path.basename(self.lmttpm), fp.write)
            print('  successfully downloaded')

        except:
            print('Error: Download failed', file=sys.stderr)
            sys.exit(1)

        #xffts
        print('Downloading a xffts file of obsid: '+str(self.obsid).zfill(6))
        try:
            with FTP(**self.config) as ftp:
                with open(self.xffts, 'wb') as fp:
                    ftp.retrbinary('RETR '+elf.localBaseDir+'xffts/'+os.path.basename(self.xffts), fp.write)
            print('  successfully downloaded')

        except:
            print('Error: Download failed', file=sys.stderr)
            sys.exit(1)

    def SearchAndDownload(self,obsid):
        self.FindRawdata(obsid)
        self.DownloadRawdata()


    ###
