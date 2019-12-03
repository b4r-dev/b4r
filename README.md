# b4rpipe
:rocket:A pipeline reduction tool for B4R/LMT data.

**The scripts are still under development. We are not responsible for the outputs now.**

**************************************************************
Python environment
**************************************************************

```terminal
$ cd b4rpipe
$ pip install pipenv
$ pipenv install
$ pipenv shell
```

You are now in the python environment.

**************************************************************
Interactive mode
**************************************************************
```terminal
$ python
$ >>> import B4Rpipe2 as Bp
$ >>> Bp.globBaseDir = '/Volumes/hdd_mac/b4r'
$ >>> Bp.globLogDir = '/Volumes/hdd_mac/b4r/logv1'
$ >>> Bp.PipelineAnalysis(87269)
```

This example shows a pipeline analysis for the data obsid 87269.

You need to specify following pathes.
#### B4Rpipe2.globBaseDir
* The path where "xffts" and "lmttpm" directorys are located.
* XFFTS binary data (e.g., xffts20181003111006.xfftsx.01) should be stored under the "xffts" directory.
* LMT antenna log data (e.g., lmttpm_2018-04-22_075858_01_0000.nc) should be stored under the "lmttpm" directory.

#### B4Rpipe2.globLogDir
* The path where outputs from B4Rpipe2 are created.
* Anywhere you like is OK.

**************************************************************
Products
**************************************************************
The script create following outputs (if possible).

* Continuum Map Qlook (Pointing offset, efficiency, etc.)
* Line (SiO) Map Qlook (Pointing offset, etc.)
* Spectrum Qlook (with auto-flag)
* MS2 (CASA readable format)
