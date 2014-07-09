from astropy import table
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import scipy.spatial as sps
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
		print key + ':\t', tree.query([x, y], distance_upper_bound = 3)
		COG[key] = np.array(cat[key]['FLUX_APER'][tree.query([x, y])[1]])
		COG_e[key] = np.array(cat[key]['FLUXERR_APER'][tree.query([x, y])[1]])
		
	#ds9.set('exit')
	return COG, COG_e
	
def SED(data, error, aperture):
	'''
	input a COG dictionary (data) with all bands (such as from obj_phot), along with a 
	maximum aperture (diameter, in arcseconds), normalize to that aperture, and return an 
	SED at that aperture
	'''
	vals = []
	errs = []
	apers = np.array([2,3,4,6,8,10,14,20,28,40,60,80,100,160])*0.06
	for key in data.keys():
		norm_val = np.interp(aperture, apers/2., data[key])
		vals.append(norm_val)
		norm_err = np.interp(aperture, apers/2., error[key])
		errs.append(norm_err)
		#plt.errorbar(bands[key], norm_val, yerr = norm_err, marker = 'x', color = 'k')
	#plt.xlabel('Wavelength (nm)')
	#plt.ylabel('Log10 of Flux')
	#plt.title('SED at aperture d = ' + str(aperture) + '"')
	#plt.yscale('log', basex = 10)
	#plt.show()
	return np.asarray(vals), np.asarray(errs)

#Generate COG dicts for a star
coords = np.array([ [2284, 2204], [3380, 2051], [2778, 1481], [2810, 3210], [3043, 1691] ])
bands = {'105': 1055.2, '125': 1248.6, '140': 1392.3, '160': 1536.9, '435': 429.7, '606': 590.7, '814': 833.3, 'K': 2200}
#check the center of the K-band
bands_wid = {'105': 265/2., '125': 284.5/2., '140': 384/2., '160': 268.3/2., '435': 103.8/2., '606': 234.2/2., '814': 251.1/2., 'K': 200.}
i = 4
'''
COG, COG_e = obj_phot(coords[i, 0], coords[i, 1])
'''
'''
#now make an SED with a bunch of different apertures for a star	
vals, errs = SED(COG, COG_e, 1.2)
for a in [.2, .4, .6, .9, 1.2]:
	vals, errs = SED(COG, COG_e, a)
	plt.errorbar(np.asarray(bands.values())/1000., vals, xerr = np.asarray(bands_wid.values())/1000., yerr = errs, marker = 'x', linestyle = 'None', label = str(a) + '"')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength (microns)')
plt.ylabel('Flux')
plt.title('SED for several apertures')
plt.legend(loc= 'best', title = 'Apertures')
plt.show()
'''

#now do the same for Y1
Y1_coords = np.array([1575, 3527])
COG, COG_e = obj_phot(Y1_coords[0], Y1_coords[1])
vals, errs = SED(COG, COG_e, 1.2)

'''
for a in [2.5, 3.5, 4.5, 5.5]:
	vals, errs = SED(COG, COG_e, a)
	plt.errorbar(np.asarray(bands.values())/1000., vals, xerr = np.asarray(bands_wid.values())/1000., yerr = errs, marker = 'x', linestyle = 'None', label = str(a) + '"')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength (microns)')
plt.ylabel('Flux')
plt.title('SED for several apertures')
plt.legend(loc= 'best', title = 'Apertures')
plt.show()
'''
#We're getting detection in all bands (!!!)

#now for the heck of it, make a COG for Y1 in all bands simultaneously
apers = np.array([2,3,4,6,8,10,14,20,28,40,60,80,100,160])*0.06
for key in COG.keys():
	plt.errorbar(apers, COG[key]/COG[key][7], marker = 'x', label = key)

plt.legend(loc = 'best')
plt.xlim([0., 8*.06])
plt.ylim([0., 1.])
plt.show()