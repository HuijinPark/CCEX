import pandas as pd
import os
import sys
from spHamilReader.spin_database.constants import * 
import spHamilReader.spin_database.SpinDict as sd

class myspindict:

	def __init__(self,s=0.,gyro=0.,q=0.,conc=.0):
		self.s = s				#: spin 		(no unit)
		self.gyro = gyro        #: gyro.. 		(rad/ms/G)
		self.q = q              #: eQ 			(barn = 10^-28 m^2)
		self.conc = conc        #: concent.. 	(%/100)
	def __repr__(self):
		try:
			base_message = f'({self.s:.1f}, {self.gyro:.4f}, {self.q:.5f}, {self.conc:.4f})'
		except TypeError:
			base_message = f'({self.s}, {self.gyro}, {self.q}, {self.conc})'
		return base_message


def isotopeData(atmsymbol,help=None):

	if help !=None:
		print("	%\n\
				%  EASY SPIN DATABASE\n\
				%\n\
				%  Column 1: #protons\n\
				%  Column 2: #nucleons\n\
				%  Column 3: radioactive *, stable -\n\
				%  Column 4: symbol\n\
				%  Column 5: name\n\
				%  Column 6: spin quantum number\n\
				%  Column 7: nuclear g factor gn\n\
				%  Column 8: natural abundance, in percent\n\
				%  Column 9: electric quadrupole moment, in barn (10^-28 m^2)\n\
				%			NaN indicates 'not measured'\n\
				%\n")

	__location__ = os.path.realpath(
		os.path.join(os.getcwd(), os.path.dirname(__file__)))
	filepath = os.path.join(__location__, 'isotopedata.txt')

	all_spins = pd.read_csv(filepath, delim_whitespace=True, header=None, comment='%',
							names=['protons', 'nucleons', 'radioactive', 'symbol', 'name', 'spin', 'g', 'conc', 'q'])
	
	stable_spins = all_spins[(all_spins['spin'] > 0) & (all_spins['conc'] > 0)]
	
	_names = stable_spins['nucleons'].astype(str) + stable_spins['symbol']
	_gyros = stable_spins['g'] / HBAR_SI * NUCLEAR_MAGNETON / 1e7 # rad/ms/G
	_quads = stable_spins['q'] # barn (10^-28 m^2) (cf. web of elements : milibarn (10^-31 m^2)
	_spins = stable_spins['spin']
	_mi = pd.MultiIndex.from_arrays([stable_spins['symbol'], _names])
	_ser = pd.Series((stable_spins['conc'] / 100).values, index=_mi)
	common_concentrations = {level: _ser.xs(level).to_dict() for level in _ser.index.levels[0]}
	common_isotopes = {level: _ser.xs(level).to_dict() for level in _ser.index.levels[0]}
	"""
	dict: Nested dict containing natural concentrations of the stable nuclear isotopes.  
	"""
	
	# Dictionary of the common isotopes. Placed in this file to avoid circular dependency
	common_isotopes = sd.SpinDict(*zip(_names, _spins, _gyros, _quads))
	"""
	SpinDict: An instance of the ``SpinDict`` dictionary, containing properties for the most of the common isotopes with 
	nonzero spin.
	The isotope is considered common if it is stable and has nonzero concentration in nature.  
	"""
	
	# electron spin
	common_isotopes['e'] = sd.SpinType('e', 1 / 2, ELECTRON_GYRO, 0)
	if atmsymbol == 'e':
		return common_isotopes['e']
	else:
		concDict = common_concentrations[atmsymbol]

		keys= list(concDict)
		isotopedata = {}

		for k in keys:
#			print(common_isotopes[k].name)
			isotopedata[k] = myspindict(*[common_isotopes[k].s, common_isotopes[k].gyro,common_isotopes[k].q,common_concentrations[atmsymbol][k]])

	return isotopedata


