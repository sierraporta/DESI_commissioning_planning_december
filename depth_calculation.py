#!/usr/bin/env python
# coding: utf-8

# ### DEPTH exposure time calculation
# 
# Based on D.Kirkby notebook

#get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import fitsio
import astropy.table
import scipy.ndimage
from astropy.table import Table


ROOT = Path('/global/cfs/cdirs/desi/spectro/redux/andes/')
assert ROOT.exists()

wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {'b': slice(0, 2751), 'r': slice(2700, 5026), 'z': slice(4900, 7781)}


## Utilities
class Spectrum(object):
    def __init__(self, stype, flux=None, ivar=None, mask=None):
        assert stype is 'full' or stype in cslice, 'invalid stype'
        self.stype = stype
        self.wave = fullwave[cslice[stype]] if stype in cslice else fullwave
        if flux is None and ivar is None:
            self._flux = np.zeros(len(self.wave))
            self.ivar = np.zeros(len(self.wave))
        elif flux is not None and ivar is not None:
            self._flux = np.asarray(flux)
            self.ivar = np.asarray(ivar)
            assert self.ivar.shape == self._flux.shape, 'flux and ivar have different shapes.'
        else:
            raise ValueError('flux and ivar must both be specified.')
        if mask is None:
            self.mask = np.zeros_like(self._flux, bool)
        else:
            self.mask = np.asarray(mask)
            assert self.mask.shape == self._flux.shape, 'flux and mask have different shapes.'
    def copy(self):
        return Spectrum(self.stype, self.flux.copy(), self.ivar.copy(), self.mask.copy())
    def __itruediv__(self, factor):
        np.divide(self.flux, factor, out=self._flux, where=factor != 0)
        self.ivar *= factor ** 2
        return self
    def __truediv__(self, factor):
        result = self.copy()
        result /= factor
        return result
    @property
    def flux(self):
        return self._flux

class CoAdd(Spectrum):
    def __init__(self, stype):
        super(CoAdd, self).__init__(stype)
        self._weighted_flux_sum = np.zeros(len(self.wave))
        self._finalized = False
    def __iadd__(self, other):
        if other.stype == self.stype:
            self_slice = slice(None, None)
        elif self.stype == 'full':
            self_slice = cslice[other.stype]
        else:
            raise ValueError(f'Cannot add "{other.stype}" to "{self.stype}".')
        self._weighted_flux_sum[self_slice] += other.ivar * other.flux
        self.ivar[self_slice] += other.ivar
        self._finalized = False
        return self
    @property
    def flux(self):
        if not self._finalized:
            np.divide(self._weighted_flux_sum, self.ivar, out=self._flux, where=self.ivar > 0)
            self._finalized = True
        return self._flux

def exposures_of_tile(tileid, night):
    paths = (ROOT / 'tiles' / str(tileid) / str(night)).glob(f'cframe-[brz][0-9]-????????.fits')
    expids = set()
    for path in sorted(paths):
        expids.add(int(path.name[-13:-5]))
    return sorted(expids)


# Read DESIMODEL throughputs for each spectrograph from a local copy obtained using:
# ```
# svn co https://desi.lbl.gov/svn/code/desimodel/tags/0.13.0/data/throughput thru13
# ```

def load_spec_thru(path='thru13'):
    thru = {}
    for camera, color in zip('brz', 'brk'):
        data = fitsio.read(f'{path}/thru-{camera}.fits', 'THROUGHPUT')
        thru[camera] = np.interp(fullwave[cslice[camera]], data['wavelength'], data['throughput'])
        plt.plot(fullwave[cslice[camera]], thru[camera], '.', c=color, ms=1)
    plt.xlabel('Wavlength [A]')
    plt.ylabel('Throughput')
    return thru

spec_thru = load_spec_thru()


# ## GFA Analysis
# 
# Read offline GFA measurements of transparency and fiber fraction from the FITS file attached to the [SV0 wiki page](https://desi.lbl.gov/trac/wiki/TargetSelectionWG/SV0):

gfa_results = fitsio.read('desi_sv0_exposures_gfa.fits', ext=1)

# ## Sky Background
# Estimate the average sky for a single exposure in phot/sec detected in each camera and incident on M1:

_sky_cache = {}

def get_sky(night, expid, specs=range(10), use_cache=True, fill_cache=True):
    """Estimate the sky spectrum for one exposure in units of phot/sec per wavelength bin.
    
    Returns a tuple (flux_inc, ivar_inc, flux_det, ivar_det) where "det" is detected phot/sec
    in each camera with flat-field corrections applied, and "inc" corrects for the average
    spectrograph throughput in each camera, then coadds over cameras.
    """
    night = str(night)
    expid = str(expid).zfill(8)
    if use_cache and (night, expid) in _sky_cache:
        return _sky_cache[(night, expid)]
    incident = CoAdd('full')
    detected = {}
    # Loop over cameras.
    for camera in 'brz':
        detected[camera] = CoAdd(camera)
        # Loop over spectrographs.
        for spec in specs:
            # Read the flat-fielded (constant) sky model in this spectrograph.
            skypath = ROOT / 'exposures' / night / expid / f'sky-{camera}{spec}-{expid}.fits'
            if not skypath.exists():
                print(f'Skipping non-existent {camera}{spec}.')
                continue
            with fitsio.FITS(str(skypath)) as hdus:
                exptime = hdus[0].read_header()['EXPTIME']
                flux = hdus['SKY'].read()
                ivar = hdus['IVAR'].read()
                mask = hdus['MASK'].read()
                # Verify that we have the expected wavelengths.
                assert np.allclose(detected[camera].wave, hdus['WAVELENGTH'].read())
                # Verify that ivar is purely statistical.
                assert np.array_equal(ivar, hdus['STATIVAR'].read())
                # Verify that the model has no masked pixels.
                assert np.all((mask == 0) & (ivar > 0))
                # Verify that the sky model is constant.
                assert np.array_equal(np.max(ivar, axis=0), np.min(ivar, axis=0))
                ##assert np.allclose(np.max(flux, axis=0), np.min(flux, axis=0))
                # There are actually small variations in flux!
                # TODO: figure out where these variations come from.
                # For now, take the median over fibers.
                detected[camera] += Spectrum(camera, np.median(flux, axis=0), np.median(ivar, axis=0))
        # Scale to the exposure time.
        detected[camera] /= exptime
        # Correct for throughput and accumulate over cameras.
        incident += detected[camera] / spec_thru[camera]
    if fill_cache:
        _sky_cache[(night, expid)] = (incident, detected)
    return incident, detected

def load(path):
    spec = {}
    with fitsio.FITS(str(path)) as hdus:
        for camera in 'brz':
            spec[camera] = hdus[camera].read()
    return spec

det_eso = load('dark_eso.fits')
det_desimodel = load('dark_desimodel.fits')

def combine_sky(night=20200314, expids=range(55355, 55368)):
    incident = CoAdd('full')
    detected = {C: CoAdd(C) for C in 'brz'}
    for expid in expids:
        exp_incident, exp_detected = get_sky(night, expid)
        incident += exp_incident
        for camera in 'brz':
            detected[camera] += exp_detected[camera]
    return incident, detected

inc_darkobs, det_darkobs = combine_sky()

def plot_dark(smoothing=125):#, save='fidsky.png'):
    plt.figure(figsize=(12, 5))
    for camera, color in zip('brz', 'brk'):
        wave = det_darkobs[camera].wave
        plt.plot(wave, scipy.ndimage.gaussian_filter1d(
            det_darkobs[camera].flux, smoothing), ':', c=color, ms=1)
        plt.plot(wave, scipy.ndimage.gaussian_filter1d(
            det_desimodel[camera], smoothing), '--', c=color, ms=1)
        plt.plot(wave, scipy.ndimage.gaussian_filter1d(
            det_eso[camera], smoothing), '-', c=color, ms=1)
    plt.plot([], [], 'k-', label='ESO')
    plt.plot([], [], 'k--', label='DESIMODEL')
    plt.plot([], [], 'k:', label='20200314')
    plt.legend()
    plt.xlabel('Wavelength [A]')
    plt.ylabel('Detected Sky [elec/sec]')
    plt.xlim(fullwave[0], fullwave[-1])
    plt.grid()
    plt.tight_layout()
        
plot_dark()

_depths = {}

def plot_tile_depth(tileid, night, darkref=det_eso, ffracref=0.56, smoothing=125):
    tileid = str(tileid)
    night = str(night)
    fig, ax = plt.subplots(2, 1, figsize=(9, 8))
    # Lookup exposures contributing to the coadd of this (tileid, night).
    expids = exposures_of_tile(tileid, night)
    # Get GFA data for these exposures.
    sel = (gfa_results['tileid'] == int(tileid)) & (gfa_results['night'] == int(night))
    gfa = gfa_results[sel]
    if not np.array_equal(expids, gfa['expid']):
        drop = set(expids) - set(gfa['expid'])
        print(f'WARNING: dropping exposures with no GFA results: {drop}.')
        expids = gfa['expid']
        assert np.all(np.diff(expids) > 0)
    nexp = len(expids)
    exptimes = np.empty((6, nexp))
    exptimes[0] = gfa['exptime']
    exptimes[1] = exptimes[0] * gfa['transparency_med'] ** 2
    exptimes[2] = exptimes[1] * (gfa['fiber_fracflux_med'] / ffracref) ** 2
    # Loop over exposures of this tile.
    for k, expid in enumerate(expids):
        inc, det = get_sky(night, expid)
        for j, camera in enumerate('brz'):
            wave = det[camera].wave
            smoothref = scipy.ndimage.gaussian_filter1d(darkref[camera], smoothing)
            if k == 0:
                ax[1].plot(wave, smoothref, 'k:', label='FIDUCIAL' if camera == 'b' else None)
            smooth = scipy.ndimage.gaussian_filter1d(det[camera].flux, smoothing)
            ax[1].plot(wave, smooth, c=f'C{k}', label=int(expid) if camera == 'b' else None)
            # Calculate the mean ratio of actual / fiducuial sky in each camera.
            #ratio = det[camera] / darkref[camera]
            #mean_ratio = np.sum(ratio.ivar * ratio.flux) / np.sum(ratio.ivar)
            mean_ratio = np.sum(smooth) / np.sum(smoothref)
            #print(k, expid, camera, 'actual/ref =', mean_ratio, mean_ratio2)
            exptimes[3 + j, k] = exptimes[2, k] / mean_ratio
        print(f'{tileid} {night} {expid} b={exptimes[3, k]:6.1f}s r={exptimes[4, k]:6.1f}s  z={exptimes[5, k]:6.1f}s')
    exptimes = np.cumsum(exptimes, axis=1)
    ax[0].plot(exptimes[0], 'c:x', label=f'Actual {exptimes[0, -1]:.0f}s')
    ax[0].plot(exptimes[1], 'c--x', label=f'Transp {exptimes[1, -1]:.0f}s')
    ax[0].plot(exptimes[2], 'c-x', label=f'FFrac {exptimes[2, -1]:.0f}s')
    for j, (camera, color) in enumerate(zip('brz', 'brk')):
        ax[0].plot(exptimes[3 + j], '-o', c=color, label=f'{camera} Sky {exptimes[3 + j, -1]:.0f}s')
    ax[0].set_xlabel('Exposure')
    ax[0].set_xticks([])
    ax[0].set_ylabel('Integrated Exposure Time [s]')
    ax[0].legend()
    ax[0].set_ylim(0, None)
    ax[0].grid()
    ax[1].legend(title=f'{tileid} {night}')
    ax[1].set_xlabel('Wavelength [A]')
    ax[1].set_ylabel('Smoothed Sky Flux [elec/sec]')
    ax[1].grid()
    plt.tight_layout()
#    plt.savefig(f'depth/depth_{tileid}_{night}.png')
    # Save results for a tabular summary.
    _depths[(tileid, night)] = exptimes

plot_tile_depth(67230, 20200314)


def determine_tile_depth2(tileid, night, expid, darkref=det_eso, ffracref=0.56, smoothing=125):
    tileid = str(tileid)
    night = str(night)
    expid = str(expid)
    # Get GFA data for these exposures.
    sel = (gfa_results['tileid'] == int(tileid)) & (gfa_results['night'] == int(night)) & (gfa_results['expid'] == int(expid))
    gfa = gfa_results[sel]
    exptimes = np.empty((6))
    exptimes[0] = gfa['exptime']
    exptimes[1] = exptimes[0] * gfa['transparency_med'] ** 2
    exptimes[2] = exptimes[1] * (gfa['fiber_fracflux_med'] / ffracref) ** 2
    inc, det = get_sky(night, expid)
    for j, camera in enumerate('brz'):
        wave = det[camera].wave
        smoothref = scipy.ndimage.gaussian_filter1d(darkref[camera], smoothing)
        smooth = scipy.ndimage.gaussian_filter1d(det[camera].flux, smoothing)
        mean_ratio = np.sum(smooth) / np.sum(smoothref)
        exptimes[3 + j] = exptimes[2] / mean_ratio
    _depths[(tileid, night)] = exptimes
    return tileid, night, expid, np.round(exptimes[3],1), np.round(exptimes[4],1), np.round(exptimes[5],1)


bdepth, rdepth,zdepth = [], [], []
for i in range(len(gfa_results)):
#    print(determine_tile_depth2(gfa_results["tileid"][i], gfa_results["night"][i], gfa_results["expid"][i]))
    bdepth.append(determine_tile_depth2(gfa_results["tileid"][i], gfa_results["night"][i], gfa_results["expid"][i])[3])
    rdepth.append(determine_tile_depth2(gfa_results["tileid"][i], gfa_results["night"][i], gfa_results["expid"][i])[4])
    zdepth.append(determine_tile_depth2(gfa_results["tileid"][i], gfa_results["night"][i], gfa_results["expid"][i])[5])


gfa_results=Table(gfa_results)
gfa_results["b_depth"]=bdepth
gfa_results["r_depth"]=rdepth
gfa_results["z_depth"]=zdepth

gfa_results.write("desi_sv0_exposures_gfa_with_depth.fits",overwrite=True)

