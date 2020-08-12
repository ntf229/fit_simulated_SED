# reads parameters from config.txt and changes .ski file
# begin line with 'SKIRT/' for skirt parameters
# use a unique parent header after the '/' to point to each parameter
# example: SKIRT/GeometricSource/scaleLength "2000 pc"

import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename")
parser.add_argument("--projectPath")
args = parser.parse_args()

tree = ET.parse(args.projectPath+'/SKIRT_files/'+args.filename)
root = tree.getroot()

# store config.txt inputs as dictionary 
d = {}
with open(args.projectPath+"/config.txt") as f:
	for line in f:
		(key, val) = line.split()
		d[key] = val

print (d)

for name, value in d.items():
	s = name.split('/')
	print('split:', s)
	print('length:', len(s))
	print(s[-1], value)

	if s[0] == 'SKIRT':

		for item in root.iter(s[1]):
			item.set(s[-1], value.replace("_", " "))


tree.write(args.projectPath+'/SKIRT_files/'+args.filename, encoding='UTF-8', xml_declaration=True)


"""


tree = ET.parse('built-ins.ski')
root = tree.getroot()

for item in root.iter('MonteCarloSimulation'):
	item.set('numPackets', "1e6")

# Wavelength grids:
for item in root.iter('sourceSystem'):
	item.set('minWavelength', "0.1 micron")
	item.set('maxWavelength', "100 micron")
	item.set('sourceBias', "0.5")
for item in root.iter('IntegratedLuminosityNormalization'):
	item.set('integratedLuminosity', '1e9 Lsun') # controls total stellar mass 
	item.set('minWavelength', '0.1 micron')
	item.set('maxWavelength', '100 micron')
for child in root.iter('GeometricSource'):
	for item in child.iter('LogWavelengthDistribution'):
		item.set('minWavelength', '0.1 micron')
		item.set('maxWavelength', '100 micron')
for child in root.iter('DustEmissionOptions'):
	for item in child.iter('LogWavelengthDistribution'):
		item.set('minWavelength', '0.1 micron')
		item.set('maxWavelength', '100 micron')
for child in root.iter('wavelengthBiasDistribution'):
	for item in child.iter('LogWavelengthDistribution'):
		item.set('minWavelength', '0.1 micron')
		item.set('maxWavelength', '100 micron')
for child in root.iter('InstrumentSystem'):
	for item in child.iter('LogWavelengthGrid'):
		item.set('minWavelength', '0.1 micron')
		item.set('maxWavelength', '100 micron')
		item.set('numWavelengths', '1000')



for child in root.iter('GeometricSource'):
	for item in child.iter('ExpDiskGeometry'):
		item.set('scaleLength', '2000 pc')
		item.set('scaleHeight', '300 pc')

for child in root.iter('GeometricMedium'):
	for item in child.iter('ExpDiskGeometry'):
		item.set('scaleLength', '2000 pc')
		item.set('scaleHeight', '300 pc')

for child in root.iter('GeometricMedium'):
	for item in child.iter('MassMaterialNormalization'):
		item.set('mass', '1e7 Msun')

for item in root.iter('FSPSSED'):
	item.set('imf', 'Chabrier')
	item.set('metallicity', '0.02')
	item.set('age', '5 Gyr')
		
for child in root.iter('InstrumentSystem'):
	for item in child.iter('FullInstrument'):
		item.set('instrumentName', 'i00')
		item.set('distance', '100 Mpc')
		item.set('inclination', '0 deg')
		item.set('azimuth', '0 deg')
		item.set('fieldOfViewX', '2e4 pc')
		item.set('numPixelsX', '250')
		item.set('centerX', '0 pc')
		item.set('fieldOfViewY', '2e4 pc')
		item.set('numPixelsY', '250')
		item.set('centerY', '0 pc')


tree.write('built-ins_2.ski', encoding='UTF-8', xml_declaration=True)

"""



