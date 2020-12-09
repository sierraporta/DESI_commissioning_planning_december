#!/usr/bin/env python

import sys                                                                                                                                            
import numpy as np
import os
import astropy.io.fits as fits
import glob

campaign = sys.argv[1] # minisv2 or sv0

outfits= './desi_'+campaign+'_exposures_gfa.fits'
expdir = '/global/cfs/cdirs/desi/spectro/redux/daily/exposures/'

if (campaign=='minisv2'):
	# assumes pid is here for all exposures
	pid   = 'b0'
	# https://desi.lbl.gov/trac/wiki/TargetSelectionWG/miniSV2
	tiles = [70002,70003,70004,70005,70006,70500,70502,70506,70508,70510,70511,70512,70513,70514]
	gfafn = '/global/cfs/cdirs/desi/users/ameisner/GFA/conditions/offline_all_guide_ccds_minisv2.fits'
	# listing all exposures (19 feb to 04 mar 2020, included)
	nights = []
	nights+= [fn.split('/')[-1] for fn in glob.glob(expdir+'202002??') if int(fn.split('/')[-1][-2:])>=19]
	nights+= [fn.split('/')[-1] for fn in glob.glob(expdir+'202003??') if int(fn.split('/')[-1][-2:])<= 4]
elif (campaign=='sv0'):
	# assumes pid is here for all exposures (taking b4 as some cframe-b{0,1,2,3} are currently missing...)
	pid = 'b4'
	# https://desi.lbl.gov/trac/wiki/TargetSelectionWG/SV0
	tiles = [65008,66000,66003,66014,66019,67142,67230,68000,68001,68002]
	gfafn = '/global/cfs/cdirs/desi/users/ameisner/GFA/conditions/offline_all_guide_ccds_SV0.fits'
	nights= ['20200314','20200315']
else:
	sys.exit('wrong campaign!')


## expid-cube_index unique identifier
gfa    = fits.open(gfafn)[1].data
gfa_eci= np.array([str(e)+'-'+str(c) for e,c in zip(gfa['expid'],gfa['cube_index'])])


allexps   = []
allnights = []
for night in nights:
	exps       = [fn.split('/')[-2] for fn in glob.glob(expdir+night+'/????????/cframe-'+pid+'-????????.fits')]
	exps       = [exp for exp in exps if fits.getheader(expdir+night+'/'+exp+'/cframe-'+pid+'-'+exp+'.fits')['tileid'] in tiles]
	exps       = [exp for exp in exps if exp!='00055587'] # 00055587 wrongly has tiled=67230
	allexps   += exps
	allnights += [night for exp in exps]
allexps   = np.array(allexps)
allnights = np.array(allnights)
print(str(len(allexps)), 'exposures reduced to cframe during '+campaign)
# list tileid (assuming b0 is here)
alltiles = np.array([fits.getheader(expdir+night+'/'+exp+'/cframe-'+pid+'-'+exp+'.fits',0)['tileid'] for night,exp in zip(allnights,allexps)])

tmp = np.argsort(allexps)
for tile,night,exp in zip(alltiles[tmp],allnights[tmp],allexps[tmp]):
	print(tile,night,exp)



# stored quantities
## cframe
cfrkeys = ['expid','night','tileid','exptime','mjdobs','tilera','tiledec']
cfrfmts = ['K',    'K',    'K',     'E',      'E',     'E',    'E']
## gfa
gfakeys,gfafmts = [],[]
for key,fmt in zip(
		['airmass','moon_sep_deg','transparency','fwhm_asec','sky_mag_ab','fiber_fracflux'],
		['E',      'E',           'E',           'E',        'E',         'E']):
	for quant in ['min','med','max']:
		gfakeys += [key+'_'+quant]
		gfafmts += [fmt]
## own
ownkeys = ['ngfa','ebv']
ownfmts = ['K','E']


# initialising
mydict = {}
for key in cfrkeys+gfakeys+ownkeys:
	mydict[key] = []


# loop on tiles
for tile in tiles:
	# exposures for that tile
	allkeep = (alltiles==int(tile))
	print(tile, allkeep.sum(), 'exposures : ',allexps[allkeep])
	# loop on exposures
	for night,exp in zip(allnights[allkeep],allexps[allkeep]):
		# exposure infos
		fn  = expdir+night+'/'+exp+'/cframe-'+pid+'-'+exp+'.fits'
		hdr = fits.getheader(fn,0)
		mydict['ebv']    += [np.nanmedian(fits.open(fn)['fibermap'].data['ebv'])]
		for key in cfrkeys:
			if   (key=='mjdobs'): mydict[key] += [hdr['mjd-obs']]
			elif (key=='night'):  mydict[key] += [int(hdr['night'])]
			else:                 mydict[key] += [hdr[key]]
		# gfa (see aaron s email 09mar2020)
		gkeep  = (np.sqrt((gfa['skyra']-hdr['tilera'])**2 + (gfa['skydec']-hdr['tiledec'])**2)<0.1)
		gkeep &= (gfa['mjd']>hdr['mjd-obs']) & (gfa['mjd']<hdr['mjd-obs']+hdr['exptime']/86400.)
		mydict['ngfa'] += [gkeep.sum()]
		gfa_exp = gfa    [gkeep]
		eci_exp = gfa_eci[gkeep]
		tmpdict = {}
		print(tile,night,exp,gkeep.sum())
		for key in gfakeys:
			gkey = key[:-4]
			quant= key[-3:]
			# first binning by expid-cube_index (see aaron s email 09mar2020)
			if (gkeep.sum()==0):
				x = -99
			else:
				x    = []
				for eic in np.unique(eci_exp):
					tmp = (eci_exp==eic)
					x  += [np.nanmedian(gfa_exp[gkey][tmp])]
			gkey = key[:-4]
			quant= key[-3:]
			# then taking min,med,max
			if   (quant=='min'): mydict[key] += [np.nanmin(x)]
			elif (quant=='med'): mydict[key] += [np.nanmedian(x)]
			elif (quant=='max'): mydict[key] += [np.nanmax(x)]
			else:                sys.exit('wrong quant!')


# sorting by increasing exposure
tmp = np.argsort(mydict['expid'])
for key in mydict.keys():
	mydict[key] = np.array(mydict[key])[tmp]

# writing fits
collist = []
for keys,fmts in zip([cfrkeys,gfakeys,ownkeys],[cfrfmts,gfafmts,ownfmts]):
	for key,fmt in zip(keys,fmts):
		print(key,fmt)
		collist.append(fits.Column(name=key,format=fmt,array=mydict[key]))
hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
hdu.writeto(outfits,overwrite=True)

# table for the wiki
hdu  = fits.open(outfits)
data = hdu[1].data
labs = ['night','expid','tileid','exptime','ebv','airmass','moon_sep_deg','transparency','fwhm_asec','sky_mag_ab','fiber_fracflux']
fmts = ['%.0f','%.0f',  '%.0f',  '%.0f',  '%.2f','%.2f',   '%.2f',        '%.2f',        '%.2f',     '%.2f',      '%.2f']
keys = [lab if lab in hdu[1].columns.names else lab+'_med' for lab in labs]
print('')
print('||'.join(['']+labs+['']))
for i in range(len(data)):
	print('||'.join(['']+[fmt%data[key][i] for key,fmt in zip(keys,fmts)]+['']))
print('')



