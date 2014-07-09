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


coords = np.array([ [2284, 2204], [3380, 2051], [2778, 1481], [2810, 3210], [3043, 1691] ])
bands = {'105': 1050, '125': 1250, '140': 1400, '160': 1600, '435': 435, '606': 606, '814': 814, 'K': 2200}
i = 0
COG, COG_e = obj_phot(coords[i, 0], coords[i, 1])

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
	
vals, errs = SED(COG, COG_e, 1.2)
for a in [.1, .2, .3, .4, .5, .6, .8, 1., 1.2, 1.4, 1.6]:
	vals, errs = SED(COG, COG_e, a)
	plt.errorbar(np.asarray(bands.values())/1000., vals, yerr = errs, marker = 'x', linestyle = 'None', label = str(a) + '"')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength (microns)')
plt.ylabel('Flux')
plt.legend(loc= 'best', title = 'Apertures')
plt.show()

print type(bands.values())

'''
#### stellar photometry, normalized to largest aperture
star_flux = c['FLUX_APER'][stars,:]
star_flux = (star_flux.T / star_flux[:,-1]).T
'''

'''
### plot curves of growth, renormalizing to 0.8"
apers = np.array([2,3,4,6,8,10,14,20,28,40,60,80,100,160])*0.06
for i in range(stars.sum()):
    yi = np.interp(0.8, apers/2., star_flux[i,:])
    #yi=1
    plt.plot(apers/2, star_flux[i,:]/yi, alpha=0.2, color='black', label = 'Star ' + str(i))
plt.title('Stellar Curves of Growth, F105')
plt.xlabel('Aperture radius (arcseconds)')
plt.xlim([0., 2.])
plt.ylabel('Fraction of max flux')
plt.ylim([0., 1.])
'''

'''
#### WFC3 IHB values
r_ihb = np.array([0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00])
f_ihb = np.array([.502, .653, .762, .813, .833, .859, .884, .897, .924, .941, .965, .975])

yi = np.interp(0.8, r_ihb, f_ihb)
plt.plot(r_ihb, f_ihb/yi, color='red', alpha=0.9, linewidth=2, zorder=-2, label = 'WFC3 IHB')
plt.legend(loc = 'best', labelspacing = 0.2, prop={'size':8})
'''
plt.show()