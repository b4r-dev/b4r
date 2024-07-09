# b4r

[![Release](https://img.shields.io/pypi/v/b4r?label=Release&color=cornflowerblue&style=flat-square)](https://pypi.org/project/b4r/)
[![Python](https://img.shields.io/pypi/pyversions/b4r?label=Python&color=cornflowerblue&style=flat-square)](https://pypi.org/project/b4r/)
[![Downloads](https://img.shields.io/pypi/dm/b4r?label=Downloads&color=cornflowerblue&style=flat-square)](https://pepy.tech/project/b4r)
[![DOI](https://img.shields.io/badge/DO/actions/workflow/status/b4r-dev/b4r/tests.yaml?label=Tests&style=flat-square)](https://github.com/b4r-dev/b4r/actions)

Reduction and analysis tools for LMT/B4R

## Installation

```shell
pip install b4r
```

## Usage

### Reduce individual data

```python
import b4r.pipe as Bp


Bp.globBaseDir = '/home/hoge/b4r'
Bp.globLogDir = '/home/hoge/b4r/logv1'
Bp.PipelineAnalysis(86420)
```

This example shows a pipeline analysis for the data obsid `86420`.

### Reduce all data (for database)

```python
import b4r.pipe as Bp

Bp.globBaseDir = '/home/hoge/b4r'
Bp.globLogDir = '/home/hoge/b4r/logv1'
Bp.PipelineAnalysisBatchRun()
```

You need to specify following pathes.

#### b4r.pipe.globBaseDir

- The path where `xffts` and `lmttpm` directories are located.
- XFFTS binary data (e.g., `xffts20181003111006.xfftsx.01`) should be stored under the `xffts` directory.
- LMT antenna log data (e.g., `lmttpm_2018-04-22_075858_01_0000.nc`) should be stored under the `lmttpm` directory.

#### b4r.pipe.globLogDir

- The path where outputs are created.
- Anywhere you like is OK.

## Products

The script create following outputs (if possible).

- Continuum Map Qlook (Pointing offset, efficiency (only for uranus), etc.)
- Line (SiO) Map Qlook (Pointing offset, etc.)
- Spectrum Qlook (with auto-flag)
- Time series spectrum of PSW data (NumPy readable format)
- GoDec calibration results (see Taniguchi et al. 2021)
- MS2 (CASA readable format)

## Correspondence

B4R 2018/2019 (obsnum<=087433) data <-> CASA MS2:

| Name | B4R | CASA MS2 |
| --- | --- | --- |
| polarization | A | (correlation or stokes) YY |
| polarization | B | (correlation or stokes) XX |
| sideband | LSB | spw 0 |
| sideband | USB | spw 1 |

## Information

- B4R webpage: http://lmtgtm.org/b4r/?lang=en
