from ftplib import FTP
import b4rpipe.specfile4 as sp
import sys
import os

class B4Rquery:

    def __init__(self,
                username,
                password,
                hosturl='ngc253.mtk.ioa.s.u-tokyo.ac.jp',
                b4rDir='/B4R',
                localBaseDir='./rawdata',
                ):

        config = {
            'host':hosturl,
            'user':username,
            'passwd':password,
            }

        self.config = config
        self.localBaseDir = localBaseDir
        self.b4rDir = b4rDir

        ftp = FTP(**config)
        self.lmttpm = ftp.nlst(b4rDir+'/lmttpm/')

        ftp = FTP(**config)
        self.xffts = ftp.nlst(b4rDir+'/xffts_links_ftp/*_xffts*')


    def FindRawdata(self,obsid):
        self.obsid = obsid

        # lmttpm
        result_lmttpm = [s for s in self.lmttpm if str(obsid).zfill(6) in s]

        if len(result_lmttpm)==1:
            print('Find one lmttpm file: ')
            print('  '+result_lmttpm[0])
            self.lmttpmFile = result_lmttpm[0]
        else:
            print('Error: No or more than one lmttpm files were found', file=sys.stderr)
            sys.exit(1)

        # xffts
        result_xffts = [s for s in self.xffts if str(obsid).zfill(6)+'_' in s]

        if len(result_xffts)>0:
            print('Find '+str(len(result_xffts))+' xffts files: ')
            for filename in result_xffts:
                print('  '+filename)

            self.xfftsFileList = result_xffts
        else:
            print('Error: No xffts files were found', file=sys.stderr)
            sys.exit(1)

    def DownloadRawdata(self):

        os.system('mkdir -p '+self.localBaseDir+'/lmttpm/')
        os.system('mkdir -p '+self.localBaseDir+'/xffts/')

        #lmttpm
        print('Downloading a lmttpm file of obsid: '+str(self.obsid).zfill(6))
        try:
            with FTP(**self.config) as ftp:
                with open(self.localBaseDir+'/lmttpm/'+os.path.basename(self.lmttpmFile), 'wb') as fp:
                    ftp.retrbinary('RETR '+self.lmttpmFile, fp.write)
            print(os.path.basename(self.lmttpmFile)+' is successfully downloaded')

        except:
            print('Error: Download failed', file=sys.stderr)
            sys.exit(1)

        #xffts
        print('Downloading xffts files of obsid: '+str(self.obsid).zfill(6))
        try:
            for xfftsFile in self.xfftsFileList:
                with FTP(**self.config) as ftp:
                    with open(self.localBaseDir+'/xffts/'+os.path.basename(xfftsFile).replace(str(self.obsid).zfill(6)+'_',''), 'wb') as fp:
                        ftp.retrbinary('RETR '+xfftsFile, fp.write)
                print(os.path.basename(xfftsFile).replace(str(self.obsid).zfill(6)+'_','')+' is successfully downloaded')

        except:
            print('Error: Download failed', file=sys.stderr)
            sys.exit(1)

    def SearchAndDownload(self,obsid):
        self.FindRawdata(obsid)
        self.DownloadRawdata()
