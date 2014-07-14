from astropy import table
from astropy.io import ascii
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
		errs.append(interp_err)
		keys.append(key)
	return np.asarray(vals), np.asarray(errs), keys
	
def cat_construct(base_cat):
	'''
	write an ASCII table with corrected photometry for all objects, in all bands 
	'''
	all_cats = {'105': table.Table.read('105_cat.fits', hdu=2), '125': table.Table.read('125_cat.fits', hdu=2), 
	'140': table.Table.read('140_cat.fits', hdu=2), '160': table.Table.read('160_cat.fits', hdu=2), 
	'435': table.Table.read('435_cat.fits', hdu=2), '606': table.Table.read('606_cat.fits', hdu=2), 
	'814': table.Table.read('814_cat.fits', hdu=2), 'K': table.Table.read('K_cat.fits', hdu=2)}

	#read base (deepest) catalog into new catalog, which provides reference object numbers used to look up objects in other catalogs
	obj_cat = table.Table(all_cats[base_cat]['NUMBER', 'X_IMAGE', 'Y_IMAGE'], names = ('NUMBER', 'X_IMAGE', 'Y_IMAGE'))
	
	print 'Building catalog...'
	for key in all_cats.keys():
		obj_cat['FLUX_APER_' + key] = all_cats[key]['FLUX_APER'] 
		# we assume that things remain in a sensible order. 
		# This will silently break if the catalogs are in a different order.
		# It can be fixed later, if absolutely necessary, by joining two tables on 'NUMBER'
		obj_cat['FLUXERR_APER_' + key] = all_cats[key]['FLUXERR_APER']
	return obj_cat
	
def cat_correct(obj_cat, max_aper, aper, corr):
	'''
	interpolate FLUX_APER_<BAND> and FLUXERR_APER_<BAND> along known apertures, 
	then given a max_aper (wherein we assume all flux to be contained), an aperture to use,
	and a star's correction dict, return corrected fluxes for all bands (forming an SED)
	'''
	import re
	apers = np.array([2,3,4,6,8,8.3333333,10,14,20,28,33.3333333,40,60,80,100,160])*0.06
	#print obj_cat
	
	#first do the easy part: interpolate FLUXERR_APER_<BAND>
	#define a function that simultaneously interpolates one point for many different COGs for a single band
	c = obj_cat.colnames
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
	
	#take correction dict and interpolate to find the value at max_aper and the value at aper
	#print corr
	corr_aper = {}
	for key in corr.keys():
		max_ap_value = np.interp(max_aper, apers, corr[key])
		ap_value = np.interp(aper, apers, corr[key])
		corr_aper[key] = max_ap_value/ap_value
		#print key, corr_aper[key]
	
	for key in corr_aper.keys():
		obj_cat['FLUX_APER_' + key] *= corr_aper[key]
		
	ascii.write(obj_cat, output = 'A2744_cat.dat')
		
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
	d_l = l[vals - errs > 0.]
	nd_l = l[vals - errs <= 0.]
	
	d_vals = vals[vals - errs > 0.]
	nd_vals = vals[vals - errs <= 0.]
	
	d_lw = lw[vals - errs > 0.]
	nd_lw = lw[vals - errs <= 0.]	

	d_errs = errs[vals - errs > 0.]
	nd_errs = errs[vals - errs <= 0.]
	
	colors = ['r', 'g', 'b', 'c', 'm', 'y']
	
	plt.scatter(l/1000., vals / correction, color = colors[i])

	#upper-limits on non-detections first
	plt.errorbar(nd_l/1000., nd_vals, xerr = nd_lw/1000., marker = 'v', color = colors[i], linestyle = 'None', markersize = 10)
	#then detections
	plt.errorbar(d_l/1000., d_vals, xerr = d_lw/1000., yerr = d_errs, marker = 'x', color = colors[i], linestyle = 'None', label = str(np.round(a, decimals = 2)) + '"')
#plt.yscale('log')

laporte = {'105': 27.5, '125': 26.32, '140': 26.26, '160': 26.25}
laporte_e = {'105': .08, '125': .04, '140': .03, '160': .04}
for i, key in enumerate(laporte.keys()):
	plt.errorbar(bands[key]/1000., 10.**((23.9 - laporte[key])/2.5), xerr = bands_wid[key]/1000., yerr = 10.**((23.9 - laporte[key] + laporte_e[key])/2.5) - 10.**((23.9 - laporte[key])/2.5) ,marker = 'o', color = 'k', linestyle = 'None')

plt.xlabel('wavelength (microns)')
plt.ylabel('flux (uJy)')
plt.title('Galaxy SED for several apertures')
plt.legend(loc = 'best', title = 'Ap. diam.')
plt.show()

#We're getting detection in all bands (!!!) (maybe)

obj_cat = cat_construct('160')
cat_correct(obj_cat, max_aper, min_aper, corr)