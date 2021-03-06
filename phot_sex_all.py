from astropy import table
from astropy.io import ascii
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import scipy.spatial as sps
from scipy.interpolate import interp1d
import pysao

def cat_compare():
	cat = {'105': table.Table.read('105_cat.fits', hdu=2), '125': table.Table.read('125_cat.fits', hdu=2), 
	'140': table.Table.read('140_cat.fits', hdu=2), '160': table.Table.read('160_cat.fits', hdu=2), 
	'435': table.Table.read('435_cat.fits', hdu=2), '606': table.Table.read('606_cat.fits', hdu=2), 
	'814': table.Table.read('814_cat.fits', hdu=2), 'K': table.Table.read('K_cat.fits', hdu=2)}

	#### select stars based on their small half-light radii

	colors = ['DarkOliveGreen', 'Indigo', 'DarkOrange', 'RosyBrown', 'Red', 'RoyalBlue', 'MediumSeaGreen', 'Yellow']
	coords = np.array([ [] ])

	for i, key in enumerate(cat.keys()):
		c = cat[key]
		stars = (c['MAG_AUTO'] < 25) & (c['FLUX_RADIUS'][:,0] < 2.5)
		print 'Number of stars in band', key, ':', sum(stars)
		stars = c[stars]
		plt.scatter(stars['Y_IMAGE'], stars['X_IMAGE'], color = colors[i], marker = '.', label = key)
		coord_i = np.vstack((stars['Y_IMAGE'], stars['X_IMAGE']))
		#print coord_i
		if i == 0:
			coords = coord_i
		else:
			coords = np.append(coords, coord_i, axis = 1)
		#if i < 2: print coords
	
		'''
		plt.scatter(c['MAG_AUTO'], c['FLUX_RADIUS'][:,0], alpha=0.5)
		plt.scatter(c['MAG_AUTO'][stars], c['FLUX_RADIUS'][stars,0], alpha=0.5, color='red')
		plt.title('Objects in SExtractor Catalog, F105')
		plt.xlabel('MAG_AUTO')
		plt.xlim([15., 35.])
		plt.ylabel('FLUX_RADIUS')
		plt.ylim([0., 250.])
		plt.show()
		'''
	plt.title('Star Positions')
	plt.legend(loc = 'best')
	plt.show()
	print np.shape(coords.T)
	coords = coords.T

#it looks like the object detection isn't as robust as I had hoped

#now I guess we just input star coords manually, get the nearest entry, and make the growth curves from that

def obj_phot(x, y, maxdist = 3):
	'''
	Given supplied coordinates (as they appear in ds9, not raw FITS format), 
	give closest object and its encircled-light COG points in all bands.
	'''
	cat = {'105': table.Table.read('105_cat.fits', hdu=2), '125': table.Table.read('125_cat.fits', hdu=2), 
	'140': table.Table.read('140_cat.fits', hdu=2), '160': table.Table.read('160_cat.fits', hdu=2), 
	'435': table.Table.read('435_cat.fits', hdu=2), '606': table.Table.read('606_cat.fits', hdu=2), 
	'814': table.Table.read('814_cat.fits', hdu=2), 'K': table.Table.read('K_cat.fits', hdu=2)}
	
	COG = {}
	COG_e = {}
	#ds9 = pysao.ds9() #can look at images, if necessary
	
	for key in cat.keys():
		#find the nearest object to the supplied coordinates in each band, and reject if it's too far away
		y_a = cat[key]['Y_IMAGE']
		x_a = cat[key]['X_IMAGE']
		tree = sps.KDTree(zip(x_a, y_a))
		#print key + ':\t', tree.query([x, y], distance_upper_bound = 3)
		COG[key] = np.array(cat[key]['FLUX_APER'][tree.query([x, y])[1]])
		COG_e[key] = np.array(cat[key]['FLUXERR_APER'][tree.query([x, y])[1]])
		
	#ds9.set('exit')
	return COG, COG_e

def uJ(data, band):
	'''
	takes a single flux measurement and a band string, and uses the ZP dict to output a
	flux in uJ (microJanskys)
	'''
	return data * 10.**(-0.4 * (ZPs[band] - 23.9 - band_corr[band]))
	
uJ_v = np.vectorize(uJ)
	
def COG_plot(COG, COG_e, max_aper):
	'''
	given COG data in many bands, make stellar COG for those bands, and return a correction
	factor per-band, from min_aper to max_aper, for all apertures in apers
	'''
	#now for the heck of it, make a COG in all bands simultaneously
	apers = np.array([2,3,4,6,8,8.3333333,10,14,20,28,33.3333333,40,60,80,100,160])*0.06
	corr = {}
	for key in COG.keys():
		max_val = np.interp(max_aper, apers, COG[key])
		max_val_err = np.interp(max_aper, apers, COG_e[key])
		corr[key] = max_val/COG[key]
		#print COG[key]/norm_val
		plt.subplot(2,1,1)
		plt.errorbar(apers, np.asarray(COG[key])/max_val, marker = style[key][1], color = style[key][0], label = key)
		plt.subplot(2,1,2)
		plt.plot(apers, np.asarray(corr[key]), marker = style[key][1], color = style[key][0], linestyle = '--')
	plt.subplot(2,1,1)
	plt.legend(loc = 'best')
	plt.xlim([0., max_aper])
	plt.ylim([0., 1.])
	plt.xlabel('Aperture Diameter (arcsec)')
	plt.ylabel('Fraction of flux')
	
	plt.subplot(2,1,2)
	plt.xlim([0., max_aper])
	plt.ylim([1., 3.])
	plt.ylabel('Correction factor to ' + str(max_aper))
	
	plt.show()
	#print corr
	return corr
	
def SED(data, error, max_aper):
	'''
	input a COG dictionary (data) with all bands (such as from obj_phot), along with a 
	maximum aperture (diameter, in arcseconds), normalize to that aperture, and return an 
	SED at that aperture
	'''
	vals = []
	errs = []
	keys = []
	apers = np.array([2,3,4,6,8,8.3333333,10,14,20,28,33.3333333,40,60,80,100,160])*0.06
	for key in data.keys():
		
		interp_val = np.interp(max_aper, apers, data[key])
		vals.append(interp_val)
		interp_err = np.interp(max_aper, apers, error[key])
		#errs.append(interp_err)
		errs.append(empty_aper_errs[key])
		keys.append(key)
	return np.asarray(vals), np.asarray(errs), keys
	
def cat_construct(base_cat):
	'''
	write an ASCII table with corrected photometry for all objects, in all bands 
	Add FLUX_AUTO in F160 and convert/correct into uJy
	Convert Counts to uJy, w/corrections
	Flux_radius 1 & 2 (PHOT_RADIUS) for F160
	'''
	all_cats = {'105': table.Table.read('105_cat.fits', hdu=2), '125': table.Table.read('125_cat.fits', hdu=2), 
	'140': table.Table.read('140_cat.fits', hdu=2), '160': table.Table.read('160_cat.fits', hdu=2), 
	'435': table.Table.read('435_cat.fits', hdu=2), '606': table.Table.read('606_cat.fits', hdu=2), 
	'814': table.Table.read('814_cat.fits', hdu=2), 'K': table.Table.read('K_cat.fits', hdu=2)}

	#read base (deepest) catalog into new catalog, which provides reference object numbers used to look up objects in other catalogs
	obj_cat = table.Table(all_cats[base_cat]['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'X_WORLD', 'Y_WORLD', 'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_RADIUS'],
	 names = ('NUMBER', 'X_IMAGE', 'Y_IMAGE', 'X_WORLD', 'Y_WORLD', 'FLUX_AUTO_160', 'FLUXERR_AUTO_160', 'FLUX_RADII'))
	 
	obj_cat['FLUX_RADIUS_1'] = obj_cat['FLUX_RADII'][:,0]
	obj_cat['FLUX_RADIUS_2'] = obj_cat['FLUX_RADII'][:,1]
	del obj_cat['FLUX_RADII']
	
	print 'Building catalog...'
	for key in all_cats.keys():
		obj_cat['FLUX_APER_' + key] = all_cats[key]['FLUX_APER']
		# **we assume that things remain in a sensible order**
		# This will SILENTLY BREAK if the catalogs are in a different order.
		# It can be fixed later, if absolutely necessary, by joining two tables on 'NUMBER'
		obj_cat['FLUXERR_APER_' + key] = all_cats[key]['FLUXERR_APER']
	
	print 'Preparing empty aperture columns...'
	obj_cat['EA_ERR_105'] = 0.
	obj_cat['EA_ERR_125'] = 0.
	obj_cat['EA_ERR_140'] = 0.
	obj_cat['EA_ERR_160'] = 0.
	obj_cat['EA_ERR_435'] = 0.
	obj_cat['EA_ERR_606'] = 0.
	obj_cat['EA_ERR_814'] = 0.
	obj_cat['EA_ERR_K'] = 0.
		
	print 'Getting errors from empty apertures...'
		
	print '\t 105...'
	wht_105 = fits.open('105_wht.fits')[0].data
	wht_105_p = np.percentile((wht_105[wht_105 != 0]).flatten(), 75.)
	print '\t 125...'
	wht_125 = fits.open('125_wht.fits')[0].data
	wht_125_p = np.percentile((wht_125[wht_125 != 0]).flatten(), 75.)
	print '\t 140...'
	wht_140 = fits.open('140_wht.fits')[0].data
	wht_140_p = np.percentile((wht_140[wht_140 != 0]).flatten(), 75.)
	print '\t 160...'
	wht_160 = fits.open('160_wht.fits')[0].data
	wht_160_p = np.percentile((wht_160[wht_160 != 0]).flatten(), 75.)
	print '\t 435...'
	wht_435 = fits.open('435_wht.fits')[0].data
	wht_435_p = np.percentile((wht_435[wht_435 != 0]).flatten(), 75.)
	print '\t 606...'
	wht_606 = fits.open('606_wht.fits')[0].data
	wht_606_p = np.percentile((wht_606[wht_606 != 0]).flatten(), 75.)
	print '\t 814...'
	wht_814 = fits.open('814_wht.fits')[0].data
	wht_814_p = np.percentile((wht_814[wht_814 != 0]).flatten(), 75.)
	print '\t K...'
	wht_K = fits.open('814_wht.fits')[0].data
	wht_K_p = np.percentile((wht_K[wht_K != 0]).flatten(), 75.)
	
	#here are the EA pixel noise errors
	EA_sig_105 = .0280
	EA_sig_125 = .0289
	EA_sig_140 = .0357
	EA_sig_160 = .0262
	EA_sig_435 = .0156
	EA_sig_606 = .0338
	EA_sig_814 = .0165
	EA_sig_K = .1107
	
	for row in obj_cat:
		#print 'Building empty aperture errors for object #' + str(row['NUMBER'])
		if row['NUMBER'] == 1:
			print '105 weight', wht_105[int(row['Y_IMAGE']), int(row['X_IMAGE'])], wht_105_p
		wht_pix_105 = wht_105[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_105'] = EA_sig_105 * np.sqrt(wht_105_p/(wht_pix_105+1))
		wht_pix_125 = wht_125[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_125'] = EA_sig_125 * np.sqrt(wht_125_p/(wht_pix_125+1))
		wht_pix_140 = wht_105[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_140'] = EA_sig_140 * np.sqrt(wht_140_p/(wht_pix_140+1))
		wht_pix_160 = wht_160[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_160'] = EA_sig_160 * np.sqrt(wht_160_p/(wht_pix_160+1))
		wht_pix_435 = wht_435[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_435'] = EA_sig_435 * np.sqrt(wht_435_p/(wht_pix_435+1))
		wht_pix_606 = wht_606[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_606'] = EA_sig_606 * np.sqrt(wht_606_p/(wht_pix_606+1))
		wht_pix_814 = wht_814[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_814'] = EA_sig_814 * np.sqrt(wht_814_p/(wht_pix_814+1))
		wht_pix_K = wht_K[int(row['Y_IMAGE']), int(row['X_IMAGE'])]
		row['EA_ERR_K'] = EA_sig_K * np.sqrt(wht_K_p/(wht_pix_K+1))

	#print obj_cat['EA_ERR_105']
	#a = raw_input('Press ENTER to continue...')

	return obj_cat
	
def cat_correct(obj_cat, max_aper, aper, corr):
	'''
	interpolate FLUX_APER_<BAND> and FLUXERR_APER_<BAND> along known apertures, 
	then given a max_aper (wherein we assume all flux to be contained), an aperture to use,
	and a star's correction dict, return corrected fluxes for all bands (forming an SED)
	'''
	import re
	apers = np.array([2,3,4,6,8,8.3333333,10,14,20,28,33.3333333,40,60,80,100,160])*0.06
	
	#first do the easy part: interpolate FLUXERR_APER_<BAND>
	#define a function that simultaneously interpolates one point for many different COGs for a single band
	c = obj_cat.colnames
	#print c
	errcols = [item for item in c if 'FLUXERR_APER_' in item]
	for col in errcols:
		f = interp1d(apers, obj_cat[col])
		del obj_cat[col]
		obj_cat[col] = f(aper)
	fluxcols = [item for item in c if 'FLUX_APER_' in item]
	#now interpolate the flux columns
	for col in fluxcols:
		f = interp1d(apers, obj_cat[col])
		del obj_cat[col]
		obj_cat[col] = f(aper)
	
	#print obj_cat['FLUX_APER_105']	
	#print corr
	
	#take correction dict and interpolate to find the value at max_aper and the value at aper
	#print corr
	corr_aper = {}
	for key in corr.keys():
		corr_val = np.interp(aper, apers, corr[key])
		fluxcolname = 'FLUX_APER_' + key
		#print fluxcolname
		#print corr_val
		corr_band = obj_cat[fluxcolname] * corr_val
		obj_cat[fluxcolname] = corr_band
		#print obj_cat[fluxcolname]
	
	#now convert fluxes and errors to uJy using uJ_v fnc (which applies extinction and ZP correction, as well)
	for errcolname in errcols:
		#print errcolname
		if 'K' in errcolname: band = 'K'
		if 'K' not in errcolname: band = errcolname[-3:]
		fluxcolname = 'FLUX_APER_' + band
		eaerrcolname = 'EA_ERR_' + band
		#print obj_cat[eaerrcolname]
		corr_flux = np.around( uJ_v(obj_cat[fluxcolname], band), 5)
		corr_fluxerr = np.around( uJ_v(obj_cat[fluxcolname] + obj_cat[errcolname], band) - uJ_v(obj_cat[fluxcolname], band), 5)
		corr_eaerr = np.around( uJ_v(obj_cat[fluxcolname] + obj_cat[eaerrcolname], band) - uJ_v(obj_cat[fluxcolname], band), 5)
		#print corr_flux
		#print corr_eaerr
		obj_cat[fluxcolname] = corr_flux
		obj_cat[errcolname] = corr_fluxerr
		obj_cat[eaerrcolname] = corr_eaerr
	
	#finally, correct the 'FLUX_AUTO_160' band
	flux_auto = uJ_v(obj_cat['FLUX_AUTO_160'], '160')
	fluxerr_auto = uJ_v(obj_cat['FLUX_AUTO_160'] + obj_cat['FLUXERR_AUTO_160'], '160') - flux_auto
	obj_cat['FLUX_AUTO_160'] = flux_auto
	obj_cat['FLUXERR_AUTO_160'] = fluxerr_auto
		
	ascii.write(obj_cat, output = 'A2744_cat_EA.dat')
	#print obj_cat
	return obj_cat
		
# =====

#Generate COG dicts for a star
coords = np.array([ [2284, 2204], [3380, 2051], [2778, 1481], [2810, 3210], [3043, 1691] ])
Y1_coords = np.array([1575, 3527])
bands = {'105': 1055.2, '125': 1248.6, '140': 1392.3, '160': 1536.9, '435': 429.7, '606': 590.7, '814': 833.3, 'K': 2200}
#check the center of the K-band
bands_wid = {'105': 265/2., '125': 284.5/2., '140': 384/2., '160': 268.3/2., '435': 103.8/2., '606': 234.2/2., '814': 251.1/2., 'K': 400./2.}
style = {'105': ['r', 'x'], '125': ['g', 'x'], '140': ['b', 'x'], '160': ['k', 'x'], '435': ['r', 'o'], '606': ['g', 'o'], '814': ['b', 'o'], 'K': ['k', 'o']}
ZPs = {'105':26.2687, '125':26.25, '140':26.46, '160':25.96, '435':25.65777, '606':26.486, '814':25.937, 'K':26.0}
band_corr = {'105':.016, '125':.012, '140':.010, '160':.008, '435':.058, '606':.040, '814':.024, 'K':.05}
empty_aper_errs = {'105':.0280, '125':.0289, '140':.0357, '160':.0262, '435':.0156, '606':.0338, '814':.0165, 'K':.1107}
apers = np.array([2,3,4,6,8,8.3333333,10,14,20,28,33.3333333,40,60,80,100,160])*0.06

#ADJUST THINGS HERE
i = 1 #which star?
max_aper = 2.0
min_aper = 0.5
# =====

#get data for a star
COG, COG_e = obj_phot(coords[i, 0], coords[i, 1])
#vals, errs, keys = SED(COG, COG_e, max_aper)
corr = COG_plot(COG, COG_e, max_aper)

'''
#now make an SED with a bunch of different apertures for a star	
for a in [.1, .2, .4, .6, .8, 1.2]:
	#print 'aper:', a
	vals, errs, keys = SED(COG, COG_e, a)
	l = np.asarray([bands[key] for key in keys])
	lw = np.asarray([bands_wid[key] for key in keys])
	#print uJ_v(vals, keys)
	#print uJ_v(vals + errs, keys) - uJ_v(vals, keys)
	#print l, lw, keys
	plt.errorbar(l/1000., uJ_v(vals, keys), xerr = lw/1000., yerr = uJ_v(vals + errs, keys) - uJ_v(vals, keys), marker = 'x', linestyle = 'None', label = str(np.round(a, decimals = 1)) + '"')
#plt.yscale('log')
plt.xlabel('wavelength (microns)')
plt.ylabel('flux (uJ)')
plt.title('Star SED for several apertures')
plt.legend(loc = 'best', title = 'ap. diam.')
plt.show()
'''

#get data for Y1
COG, COG_e = obj_phot(Y1_coords[0], Y1_coords[1])
#COG_plot(COG, COG_e, max_aper) #galaxy COG isn't very instructive

#now do the same for Y1
for i, a in enumerate([0.5]):
	print 'aper:', a
	vals, errs, keys = SED(COG, COG_e, a)
	
	#now correct the vals by a factor according to the apertures
	#correction = np.asarray([np.interp(a, apers, corr[key]) for key in keys])
	correction = []
	for key in keys:
		#print 'band:', key
		corr_i = np.interp(a, apers, corr[key])
		#print 'correction:', corr_i
		correction.append(corr_i)
	for row in zip(keys, correction): print row
	
	l = np.asarray([bands[key] for key in keys])
	lw = np.asarray([bands_wid[key] for key in keys])
	vals = uJ_v(vals, keys) * correction
	errs = uJ_v(vals + errs, keys) - uJ_v(vals, keys)
	
	#prepare to plot the detections and non-detections separately
	d_l = l[vals > 0.]
	nd_l = l[vals <= 0.]
	
	d_vals = vals[vals > 0.]
	nd_vals = vals[vals <= 0.]
	
	d_lw = lw[vals > 0.]
	nd_lw = lw[vals <= 0.]	

	d_errs = errs[vals > 0.]
	nd_errs = errs[vals <= 0.]
	
	colors = ['r', 'g', 'b', 'c', 'm', 'y']
	
	#plt.scatter(l/1000., vals / correction, color = colors[i])

	#upper-limits on non-detections first
	plt.errorbar(nd_l/1000., nd_vals, xerr = nd_lw/1000., marker = 'v', color = colors[i], linestyle = 'None')
	#then detections
	plt.errorbar(d_l/1000., d_vals, xerr = d_lw/1000., yerr = d_errs, marker = 'x', color = colors[i], linestyle = 'None', label = str(np.round(a, decimals = 2)) + '"')
	#now plot data from IRAC channel (data is already in uJy)
	plt.errorbar([3.55, 4.493], [.0961, .3770], xerr = [.75/2, 1.015/2], yerr = [.0769, .1086], marker = 'x', color = colors[i], linestyle = 'None')
#plt.yscale('log')

laporte = {'105': 27.5, '125': 26.32, '140': 26.26, '160': 26.25}
laporte_e = {'105': .08, '125': .04, '140': .03, '160': .04}
for i, key in enumerate(laporte.keys()):
	plt.errorbar(bands[key]/1000., 10.**((23.9 - laporte[key])/2.5), xerr = bands_wid[key]/1000., yerr = 10.**((23.9 - laporte[key] + laporte_e[key])/2.5) - 10.**((23.9 - laporte[key])/2.5) ,marker = 'o', color = 'k', linestyle = 'None', label = 'Laporte' if i == 0 else '')
#now plot Laporte's IRAC measurements: first detection, then non-detection
plt.errorbar(4.493, 10**((23.9 - 25.16)/2.5), xerr = 1.015/2, yerr = 10**((23.9-25.16-.16)/2.5) - 10**((23.9 - 25.16)/2.5), marker = 'x', color = 'k')
plt.errorbar(3.55, 10**((23.9 - 25.48)/2.5), xerr = .75/2, marker = 'v', color = 'k')
		
plt.xlabel('wavelength (microns)')
plt.ylabel('flux (uJy)')
plt.title('Galaxy SED for several apertures')
plt.legend(loc = 'best', title = 'Ap. diam.')
#plt.yscale('log')
plt.show()

print 'Final photometry in uJy for Abell 2744-Y1'
for row in zip(keys, vals, errs): print row[0] + '\t' + str(np.round(row[1], decimals = 4)) + '\t +/- \t' + str(np.round(row[2], decimals = 4)) + ' uJy'
print 'IRAC1' + '\t' + str(.0961) + '\t +/- \t' + str(.0769) + ' uJy'
print 'IRAC1' + '\t' + str(.3770) + '\t +/- \t' + str(.1086) + ' uJy'

#We're getting detection in all bands (!!!) (maybe)

obj_cat = cat_construct('160')
cat_correct(obj_cat, max_aper, min_aper, corr)