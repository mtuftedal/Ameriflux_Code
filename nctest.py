import os
from netCDF4 import Dataset
import numpy as np

# Path to the data stored on my local computer. Combines data to the path to be found.
FLUX_path=(r'/home/mtuftedal/Documents/ARM-to-Ameriflux/sgpfluxE39.b1/')
FLUX_files=os.listdir(FLUX_path)
print (FLUX_files)

# Sorts the files in numerical order. Splits the name up to verify a .nc
# file exists. 
for files in sorted(FLUX_files):
	name=files.split('.')
	date=name[4]
	if '.nc' in files:
		fpath_FLUX=os.path.join(FLUX_path,files)
		FLUX_data=Dataset(fpath_FLUX,mode='r')
		time=FLUX_data['time'][0:]
		time_bounds=FLUX_data['time_bounds_ecor'][0:]

		# ECORSF Data. Ameriflux data variable assigned for consistency. 
		FC=FLUX_data['co2_flux'][0:]
		# FCH4?
		#h20_flux is unsure of Ameriflux Variable
		CO2=FLUX_data['co2_molar_fraction'][0:]
		CO2_MIXING_RATIO=FLUX_data['co2_mixing_ratio'][0:]
		#H2O=FLUX_data['h2o_molar_fraction'][0:]
		H2O_MIXING_RATIO=FLUX_data['h2o_mixing_ratio'][0:]
		#CH4=FLUX_data['ch4_molar_fraction'][0:]
		#CH4_MIXING_RATIO=FLUX_data['ch4_mixing_ratio'][0:]
		TAU=FLUX_data['momentum_flux'][0:]

		



		vwc_1=FLUX_data['vwc_1'][0:]
		WS=FLUX_data['mean_wind'][0:]

		print (FLUX_data)
		#print(np.shape(time),np.shape(FC))
		print ((time_bounds))


