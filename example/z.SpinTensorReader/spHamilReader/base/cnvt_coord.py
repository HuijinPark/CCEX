from spHamilReader.base.pyform import _error
import numpy as np

def coordinate_converter( opt, unit_cell_parameter, unit_cell_information, latt_constant = 1.000000):

	unit_cell_parameter = np.asarray(unit_cell_parameter) * latt_constant
	unit_cell_information = unit_cell_information
	
	if type(unit_cell_information[0][0]) == str:pass;
	else:
		_error('The type of unit cell information is wrong')

	converted_information = []

	for i, atmdat in enumerate(unit_cell_information):
		print(atmdat)
		preconvert = np.asarray([atmdat[1],atmdat[2],atmdat[3]])
		if opt == 'frac2cart':
			postconvert = frac2cart(unit_cell_parameter,preconvert)
		elif opt == 'cart2frac':
			postconvert = cart2frac()
		else:
			_error('opt is wrong')

		dat = [atmdat[0],postconvert[0],postconvert[1],postconvert[2]]
		converted_information.append(dat)

	return converted_information

def frac2cart(unit_cell_parameter, frac_atmdat):
	return (frac_atmdat.reshape(1,3) @ unit_cell_parameter).reshape(3)

def cart2frac(unit_cell_parameter, cart_atmdat):
	frac_atmdat =  cart_atmdat.reshape(1,3) @ np.linalg.inv(unit_cell_parameter)
	if frac_atmdat[0] <= 1.000000 and frac_atmdat[1] <= 1.000000 and frac_atmdat[2] <= 1.000000:
		return frac_atmdat.reshape(3)
	else: _error('During converting cart2frac, we found that fractional coord is larger than 1.000')

