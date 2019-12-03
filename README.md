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

You need to specify following pathes.
* B4Rpipe2.globBaseDir
The path where "xffts" and "lmttpm" directorys are located.
