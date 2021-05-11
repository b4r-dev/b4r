# b4rpipe

[![](https://img.shields.io/pypi/v/b4rpipe.svg?label=PyPI&style=flat-square)](https://pypi.org/pypi/b4rpipe/)
[![](https://img.shields.io/pypi/pyversions/b4rpipe.svg?label=Python&color=yellow&style=flat-square)](https://pypi.org/pypi/b4rpipe/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?label=License&style=flat-square)](LICENSE)

A pipeline reduction tool for B4R/LMT data.

**The scripts are still under development. We are not responsible for the outputs now.**

**************************************************************
Installation
**************************************************************

```terminal
$ pip install b4rpipe
```

**************************************************************
Usage
**************************************************************

### Reduce individual data

```terminal
$ python
$ >>> import b4rpipe as Bp
$ >>> Bp.globBaseDir = '/home/hoge/b4r'
$ >>> Bp.globLogDir = '/home/hoge/b4r/logv1'
$ >>> Bp.PipelineAnalysis(86420)
```

This example shows a pipeline analysis for the data obsid 86420.

### Reduce all data (for database)

```terminal
$ python
$ >>> import b4rpipe as Bp
$ >>> Bp.globBaseDir = '/home/hoge/b4r'
$ >>> Bp.globLogDir = '/home/hoge/b4r/logv1'
$ >>> Bp.PipelineAnalysisBatchRun()
```

You need to specify following pathes.
#### b4pipe.globBaseDir
* The path where "xffts" and "lmttpm" directorys are located.
* XFFTS binary data (e.g., xffts20181003111006.xfftsx.01) should be stored under the "xffts" directory.
* LMT antenna log data (e.g., lmttpm_2018-04-22_075858_01_0000.nc) should be stored under the "lmttpm" directory.

#### b4rpipe.globLogDir
* The path where outputs are created.
* Anywhere you like is OK.

**************************************************************
Products
**************************************************************
The script create following outputs (if possible).

* Continuum Map Qlook (Pointing offset, efficiency, etc.)
* Line (SiO) Map Qlook (Pointing offset, etc.)
* Spectrum Qlook (with auto-flag)
* MS2 (CASA readable format)

**************************************************************
Data query and download (only for internal use now)
**************************************************************
If you are in the NAOJ or IoA (U. Tokyo) local network, you can access the B4R ftp server.

```terminal
$ python
$ >>> import b4rpipe as Bp
$ >>> Bp.PipelineAnalysis(86420,DataDownload=True,username='hogehoge',password='*****')
```

Then "raw" and "calibrated" directory appears at the current directory.

```terminal
$ ls
raw calibrated
$ ls raw
lmttpm xffts
$ls calibrated
86420
```

**************************************************************
Correspondence
**************************************************************
+ B4R 2018/2019 (obsnum<087433) data <-> CASA MS2
| Name | B4R | CASA MS2 |
| --- | --- | --- |
| polarization | A | YY |
| polarization | B | XX |
| sideband | LSB | 0 |
| sideband | USB | 1 |

**************************************************************
Information
**************************************************************
* B4R webpage: http://lmtgtm.org/b4r/?lang=en
* Contact: Yuki Yoshimura
  (email: astro.yoshimura(_at_)gmail.com)
