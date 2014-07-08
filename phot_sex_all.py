from astropy import table
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

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

#now try clustering with DBSCAN
db = DBSCAN(eps = 2., min_samples = 8).fit(coords)
core_samples = db.core_sample_indices_
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

import pylab as pl

unique_labels = set(labels)
colors = pl.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
	if k == -1:
		col = 'k'
		markersize = 6
	class_members = [index[0] for index in np.argwhere(labels == k)]
	cluster_core_samples = [index for index in core_samples if labels[index] == k]
	for index in class_members:
		x = coords[index]
		if index in core_samples and k != -1:
			markersize = 14
		else:
			markersize = 6
		plt.plot(x[0], x[1], 'o', markerfacecolor = col, markeredgecolor = 'k', markersize = markersize)
		
plt.show()

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