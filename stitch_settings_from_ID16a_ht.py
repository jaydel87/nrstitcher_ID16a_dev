#!/data/id16a/inhouse1/sware/NRstitcher/pyenv/bin/python3

"""
stitch_settings_from_ID16a_ht.py:

This script is dedicated to holotomographic data acquired on the beamline ID16A at the ESRF, and is used to generate a stitch settings file to be used as input to nr_stitcher.py. 
The stitch settings file is generated based on the template default_stitch_settings.py provided with nr_stitcher, with parameters mofdified based on user-selected options.
After adding the NRStitcher directory to the path, this script should be run from the directory where the generated stitch settings file is to be saved.

"""

import sys
import argparse
import glob
import re
import os
import h5py
import numpy as np
import shutil
from logFileParser import logFileParser
from itertools import compress

from pi2py2 import *
from default_stitch_settings import *

__author__ = "Jayde Livingstone"
__credits__ = ["Jayde Livingstone", "Arttu Miettenen", "Yervand Azatian"]
__maintainer__ = "Jayde Livingstone"
__email__ = "jayde.livingstone@gmail.com"


argparser = argparse.ArgumentParser(description='Creates input file for nr_stitcher.py program from id16a holotomography scan settings.')
argparser.add_argument('--base_dir', type=str, help='Required. Base directory containing sub-directories for individual scans and reconstructed volumes.')
argparser.add_argument('--vol_dir', type=str, default="volfloat/", help='Required. Sub-directory containing reconstructed volumes.')
argparser.add_argument('--volumes', nargs='*', type=str, action='store', default='*.raw', help='Required. Names (separated by a space) of reconstructed volumes to be stitched, or, if there are many volumes a single input with regular expressions can be used e.g., sample_*.raw.')
argparser.add_argument('--distances', nargs='*', type=int, action='store', default=None, help='Optional. Integers (separated by a space) corresponding to distances used in holotomographic acquisition. Default is [1 2 3 4].')
argparser.add_argument('--bitdepth', default=None, type=str, help='Optional. Bit depth of output .raw volumes if conversion is required. Default is none, no conversion will be performed. Options are "pyhst" (32-bit float), "8bits", "16bits" and "32bits" (int).')
argparser.add_argument('--voxelsize', default=None, type=str, help='Optional. Give desired voxel size (in nm) if volumes are to be rescaled (same value applies to all volumes).')
argparser.add_argument('--vrange', default=None, nargs=2, type=int, help='Optional. Minimum and maximum greyscale values. If not given they will be automatically calculated. Values are applied to all volumes')
argparser.add_argument('-b', '--binning', default=1, type=int, help='Optional. Binning that is to be applied before stitching. Recommended for optimising the stitch settings. Default is 1 (no binning).')
argparser.add_argument('-n', '--name', default='stitched', type=str, help='Sample name. Default is "stitched".')
argparser.add_argument('-o', '--output', default='stitch_settings.txt', type=str, help='Optional. Name of stitch_settings file. Default is "stitch_settings.txt". Will be saved in the current working directory.')

args = argparser.parse_args()
if args.distances == None:
	args.distances = [1, 2, 3, 4]
args.distances.sort()

	
def get_h5file(volume_name):
	"""
	Gets the name of the h5 file associated with the experiment base directory.
	
	Arguments:
	None
	
	Returns:
	Name of the h5 file as a string.
	"""
	
	h5path = args.base_dir + "../"
	h5_filename = glob.glob(h5path + "*.h5")[0]
	
	return h5_filename
	

class ht_vol:
	def __init__(self, args, volName):
		
		self.pi = Pi2()
		
		dt = {"float32" : "float32", "32bits" : "uint32", "16bits" : "uint16", "8bits" : "uint8"}
				
		self.volName = os.path.join(args.base_dir, args.vol_dir, volName)
		infoFile = os.path.join(args.base_dir, args.vol_dir, (volName + ".info"))
		self.info = self.get_vol_info(infoFile)
		self.slices = int(self.info["NUM_Z"])
		self.rows = int(self.info["NUM_Y"])
		self.cols = int(self.info["NUM_X"])
		self.voxelSize = float(self.info["voxelSize"])
		self.vmin = float(self.info["ValMin"])
		self.vmax = float(self.info["ValMax"])
		self.byteorder = self.info["BYTEORDER"].strip()
		
		self.dimensions = [self.cols, self.rows, self.slices]
		
		self.fromType = "float32"
		if args.bitdepth == None:
			args.bitdepth = self.fromType
			
		self.toType = dt[args.bitdepth]
		
		if args.vrange is not None:
			self.tomin = args.vrange[0]
			self.tomax = args.vrange[1]
		else:
			self.tomin = self.vmin
			self.tomax = self.vmax
					
		self.rescaleTo = args.voxelsize
		if self.rescaleTo != None:
			self.zoom = self.voxelSize / (float(self.rescaleTo) / 1000.)
		else:
			self.zoom = 1
		
		strbits = args.bitdepth
		
		self.modified = False
		
		self.x = round(self.cols * self.zoom)
		self.y = round(self.rows * self.zoom)
		self.z = round(self.slices * self.zoom)
		
		self.shape = str(self.y) + "x" + str(self.x) + "x" + str(self.z)
		
		self.savePath = os.path.join(args.base_dir, "volraw/")
		self.saveFile = self.volName.replace(".vol", strbits).split("/")[-1]
		if self.zoom != 1:
			self.saveFile += "_rescaled_"+str(self.rescaleTo) + "nm"
		self.saveFileName = self.saveFile + "_" + self.shape + ".raw"


	def get_vol_info(self, infoFile):
		"""
		Get the information contained with the .vol.info file.
		
		Arguments:
		infoFile: The name of the info file as a string.
		
		Returns:
		A dictionary with the property name as the KEY and its value as the VALUE.
		"""
		
		info = {}
		with open(infoFile) as f:
			for line in f.readlines():
				try:
					key, value = line.replace(" ", "").replace("\n", "").split("=")
					info[key] = value
				except:
					continue
				
		return info
	
	
	def check_raw(self):
		"""
		Check if a processed image corresponding to selected options already exists or not.
		
		Arguments:
		None
		
		Returns:
		True or False
		"""
		
		infoFileName = os.path.join(self.savePath, self.saveFileName) + ".info"
		if os.path.isfile(os.path.join(self.savePath, self.saveFileName)) and os.path.isfile(infoFileName):
			raw_info = self.get_vol_info(infoFileName)
			
			raw_vmin = int(float(raw_info["ValMin"]))
			raw_vmax = int(float(raw_info["ValMax"]))
			
			if raw_vmin == int(self.tomin) and raw_vmax == int(self.tomax):
				print(os.path.join(self.savePath, self.saveFileName) + " already exists. Skipping.")
				return True
				
		return False
			
			
	def read_vol(self):
		"""
		Opens a .vol volume as a pi2 image.
		
		Arguments:
		None
		
		Returns:
		
		"""
		
		print("Reading volume " + self.volName)
		self.volfloat = self.pi.read(self.volName)
		
		
	def create_vol(self, dtype):
		"""
		Create a new pi2 image of specified datatype and dimensions.
		
		Arguments:
		dtype: Desired datatype.
		
		Returns:
		vol: Empty pi2 image.
		"""
		
		if dtype == "float32":
			vol = self.pi.newimage(ImageDataType.FLOAT32, self.y, self.x, self.z)
		elif dtype == "uint32":
			vol = self.pi.newimage(ImageDataType.UINT32, self.y, self.x, self.z)
		elif dtype == "uint16":
			vol = self.pi.newimage(ImageDataType.UINT16, self.y, self.x, self.z)
		elif dtype == "uint8":
			vol = self.pi.newimage(ImageDataType.UINT8, self.y, self.x, self.z)
					
		return vol
		
	
	def bit_conversion(self):
		"""
		Convert volume from original datatype to new datatype with specified range.
		
		Arguments:
		None
		
		Returns:
		
		"""
		
		intbits = {"uint32" : 32, "uint16" : 16, "uint8" : 8}
		print("Converting to " + str(intbits[self.toType]) + " bits.")

		self.vmax = self.tomax
		self.vmin = self.tomin
		self.vmin -= 1 / (2 ** intbits[self.toType] - 1) # Zero-valued pixels (eg. pixels outside the reconstructed cylinder) should retain their zero-value after the conversion
		
		vrange = self.vmax - self.vmin
		
		self.pi.replace(self.volfloat, 0, self.vmin)
		self.pi.subtract(self.volfloat, self.vmin)
		self.pi.divide(self.volfloat, vrange)
		self.pi.max(self.volfloat, 0)
		self.pi.min(self.volfloat, 1)
		self.pi.multiply(self.volfloat, (2 ** intbits[self.toType] - 1))
		
		if self.toType == "uint32":
			self.pi.convert(self.volfloat, ImageDataType.UINT32)
		elif self.toType == "uint16":
			self.pi.convert(self.volfloat, ImageDataType.UINT16)
		elif self.toType == "uint8":
			self.pi.convert(self.volfloat, ImageDataType.UINT8)
		
		self.modified = True
		

	def rescale(self):
		"""
		Rescale (modify the dimensions of) the volume.
		
		Arguments:
		None
		
		Returns:
		
		"""
		
		print("Rescaling to " + str(self.rescaleTo) + " nm pixel/voxel size.")
		tmp = self.create_vol(self.toType)
		self.pi.scale(self.volfloat, tmp, [0,0,0], False, "Nearest")
		self.volfloat.set_data(tmp.get_data())
		del tmp
		self.modified = True
		
	
	def reslice(self, direction):
		"""
		Reslice the volume.
		
		Arguments:
		direction: The face in the cube which will become the first slice in the output volume. Options are Top, Bottom, Left or Right.
		
		Returns:
		
		"""
		
		tmp = self.pi.newlike(self.volraw)
		pi.reslice(self.volraw, tmp, direction)
		self.volraw = tmp.get_data()
		self.saveFile += "_resliced_"+direction
		self.modified = True
				
		
	def write_to_raw(self):
		"""
		Save the processed volume in the raw format.
		
		Arguments:
		None
		
		Returns:
		
		"""
		
		print("Writing to " + os.path.join(self.savePath, self.saveFileName))
		
		self.pi.writeraw(self.volfloat, os.path.join(self.savePath, self.saveFile))			
		self.write_raw_info()
		
						
	def write_raw_info(self):
		"""
		Write a .raw.info file containing the parameters of the processed volume.
		
		Arguments:
		None
		
		Returns:
		
		"""
		
		dims = self.volfloat.get_dimensions()
		infoFile = self.savePath + self.saveFile + "_" + str(dims[0]) + "x" + str(dims[1]) + "x" +  str(dims[2]) + ".raw.info"
		
		with open(infoFile, 'w') as f: 
			f.write("NUM_X = " + str(dims[0]) + "\n")
			f.write("NUM_Y = " + str(dims[1]) + "\n")
			f.write("NUM_Z = " + str(dims[2]) + "\n")
			f.write("voxelSize = " + str(self.voxelSize / self.zoom) + "\n")
			f.write("BYTEORDER = " + str(self.byteorder) + "\n")
			f.write("ValMin = " + str(self.tomin) + "\n")
			f.write("ValMax = " + str(self.tomax) + "\n")
	
			
	def process_vol(self):
		"""
		Calls the functions required to process the volume prior to stitching.
		
		Arguments:
		None
		
		Returns:
		The full path to the processed volume.
		"""
		
		raw_exists = self.check_raw()
		
		if not raw_exists:
			self.read_vol()
		
			if self.toType != "float32":
				self.bit_conversion()
			
			if self.zoom != 1:
				self.rescale()
		
			if self.modified:	
				self.write_to_raw()
			else:
				return self.volName	
				
		return os.path.join(self.savePath + self.saveFileName)

class h5:
	
	def __init__(self, h5_filename, volume_name):
		self.filename = h5_filename
		self.file = h5py.File(self.filename, "r")
		self.keys = list(self.file.keys())
		
		first_dist = args.distances[0]
		scan_suffix = f"_{first_dist}_"
		
		mask = [x.endswith(scan_suffix) for x in self.keys]
		entry_list = list(compress(self.keys, mask))
		
		found = False

		for entry in entry_list:
			prefix = entry.split(' ')[-1].replace(scan_suffix, '')
			if prefix in volume_name:
				self.root = self.file[entry]
				print(volume_name)
				found = True

		if found == False:
			for entry in entry_list:
				tmp = entry.split(' ')[-1]
				prefix = ''.join(tmp.rsplit(scan_suffix, 1))

				if prefix in volume_name:
					self.root = self.file[entry]
					print(volume_name)

		
	def toList(self, group):
		"""
		Write h5 file entries to a list.
		
		Arguments:
		group: Required group in h5 file
		
		Returns:
		items: names of datasets
		data: values in datasets
		"""
		
		# Names
		items = []
		for item in group.items():
			items.append(str(item[0]))
			
		data = []
		for item in group.values():
			# Data is saved as a string with spaces between values. Split on spaces to get the values as strings.
			raw_data = item[()].astype(str).split(' ')
			# Remove null strings
			filtered_data = list(filter(None,raw_data))
			data.append(filtered_data)
			
		return items, data
			
	def get_motor_positions(self):
		"""
		Get the parameters saved under sample/positioners in the experiment .h5 file.
		
		Arguments:
		None
		
		Returns:
		motor_positions_dict: Dictionary containing the positioner names as KEYS and their positions as VALUES 
		"""
		
		motor_info = self.root['sample']['positioners']
		h5_items, h5_data = self.toList(motor_info)
		names = h5_data[0] # List of motor names as strings
		values = [float(x) for x in h5_data[1]] # List of motor positions as floats
		# Create dictionary of motor positions
		motor_positions_dict = {}
		for idx,key in enumerate(names):
			motor_positions_dict[key] = values[idx]
			
		return motor_positions_dict
		
		
	def get_tomo_parameters(self):
		"""
		Get the parameters saved under the heading "TOMO" in the experiment .h5 file.
		
		Arguments:
		None
		
		Returns:
		tomo_parameters_dict: Dictionary containing the parameter names as KEYS and their values as VALUES
		"""
		
		h5_items, h5_data = self.toList(self.root['TOMO'])
		tomo_parameters_dict = {}
		for idx,key in enumerate(h5_items):
			if key != "FTOMO_PAR":
				tomo_parameters_dict[key] = float(h5_data[idx][0])
		
		return tomo_parameters_dict
	
	
def get_sax_say(su, sv):
	"""
	Calculate the sax and say values from the su and sv positions.
	
	Arguments:
	su: Position of motor su, read from experiment .h5 file
	sv: Position of motor sv, read from experiment .h5 file
	
	Returns:
	sax: Position of virtual motor sax
	say: Position of virtual motor say
	"""
	
	a = 21.5 * np.pi/180.

	sax = -np.sin(a) * su + np.cos(a) * sv
	say = np.cos(a) * su + np.sin(a) * sv
	
	return sax, say
	
	
def get_delta_sx(sx, sx0, zoom):
	"""
	Find the shift in sx which would be required to magnify (or demagnify) an image by a given factor.
	This is particularly useful for finding sx and corresponding sy and sz values for a rescaled image.
	the equation has been deduced from ...
	
	Arguments:
	sx: The position of motor sx read from the h5 file for the original image
	sx0: The fixed sx0 value (depends on beam energy), read from the experiment .h5 file 
	zoom: The magnification factor
	
	Returns:
	delta_sx: The shift in sx to be applied
	"""
	
	sx_ = (sx - sx0) * (1. / zoom) + sx0
	
	delta_sx = sx_ - sx
	
	return delta_sx
	
	
def get_delta_sy(energy, delta_sx):
	"""
	Find the shift in sy which is required to correct for the slope of the beam when moving sx.
	The angle values were obtained from the beamline macros. 
	Check with Peter Cloetens for eventual updates to these values.
	The values were last updated in the beamline macros on 31/10/2020 and were up to date as of 21/06/2023.
	Currently not used in _main_, but has been written in case it is required in future updates.
	
	Arguments:
	energy: Beam energy in keV
	delta_sx: The shift, or translation to be performed with sx
	
	Returns:
	delta_sy: shift in sy which corresponds to the given shift in sx
	"""
	
	if int(energy) == 17:
		angle = 28.866
	else:
		angle = 15.871
		
	delta_sy = -delta_sx * angle / 1000
	
	return delta_sy
	
def get_delta_sz(energy, delta_sx):
	"""
	Find the shift in sz which is required to correct for the slope of the beam when moving sx.
	The angle values were obtained from the beamline macros. 
	Check with Peter Cloetens for eventual updates to these values.
	The values were last updated in the beamline macros on 31/10/2020 and were up to date as of 21/06/2023.
	
	Arguments:
	energy: Beam energy in keV
	delta_sx: The shift, or translation to be performed with sx
	
	Returns:
	delta_sz: shift in sz which corresponds to the given shift in sx
	"""
	
	if int(energy) == 17:
		angle = 29.389
	else:
		angle = 15.415
		
	delta_sz = -delta_sx * angle / 1000
	
	return delta_sz
	
	
def motorpos_to_pixels(value, voxel_size):	
	# Converts motor position in mm to pixels
	return round(value / voxel_size)
	
	
def get_vol_info(infoFile):
	"""
	Extracts information from .info file.
	
	Arguments:
	infoFile: The name and path of the .info files to be read
	
	Returns:
	info: Dictionary containing the contents of the file, 
	where the left column, or parameter name is the KEY and the right column is the VALUE
	"""
	
	info = {}
	with open(infoFile) as f:
		for line in f.readlines():
			try:
				key, value = line.replace(" ", "").split("=")
				info[key] = value
			except:
				continue
			
	return info
	

def main():
	
	# Get volumes from command line input
	if len(args.volumes) > 1:
		# Do this if entire name of individual volumes given
		volNames = args.volumes
		floatVols = [os.path.join(args.base_dir, args.vol_dir) + x for x in volNames]
	else:
		# Do this if script is to take all volumes containing the given strings, can contain wildcards
		floatVols = sorted(glob.glob(os.path.join(args.base_dir, args.vol_dir) + args.volumes[0]))
		print(floatVols)
		volNames = [os.path.basename(x) for x in floatVols]
		
	rawVols = []
	
	for volName in volNames:
		# Process the float volumes (convert to raw if not already done)
		vol = ht_vol(args, volName)
		rawVol = vol.process_vol()
		rawVols.append(rawVol)
		
	
	floatInfoFiles = [vol + ".info" for vol in floatVols]	
	rawInfoFiles = [vol + ".info" for vol in rawVols]
	
	line_to_write = ""
	
	for i,vol in enumerate(rawVols):
		
		# Read the .info files for the .raw and original .vol volumes
		rawInfo = get_vol_info(rawInfoFiles[i])
		floatInfo = get_vol_info(floatInfoFiles[i])
		
		# Dimensions in number of pixels of the volumes
		Xdim = float(rawInfo['NUM_X'])
		Ydim = float(rawInfo['NUM_Y'])
		Zdim = float(rawInfo['NUM_Z'])
		
		# Assume standard image size pf 3216 x 3216 x 3216 px^3
		Xdim_ref = 3216
		Ydim_ref = 3216
		Zdim_ref = 3216
	
		# Get scan settings from .h5 file (the one from the experiment, containing all the scans)
		settings_file = get_h5file(vol)
		scan_settings = h5(settings_file, vol)
		motor_positions = scan_settings.get_motor_positions()
		tomo_params = scan_settings.get_tomo_parameters()
	
		sx = motor_positions['sx']
		sy = motor_positions['sy']
		sz = motor_positions['sz']
		su = motor_positions['su']
		sv = motor_positions['sv']
		sax, say = get_sax_say(su, sv)
		
		spy = motor_positions['spy'] / 1000. # spy and spz stored in um, convert to mm
		spz = motor_positions['spz'] / 1000. # spy and spz stored in um, convert to mm
		sx0 = tomo_params['sx0']
		
		orig_vox_size = float(floatInfo['voxelSize']) / 1000. # Voxel size stored in um, convert to mm 
		vox_size = float(rawInfo['voxelSize']) / 1000. # Voxel size stored in um, convert to mm
		
		Z = sz
		Y = -sax
		X = -say
			
		# Magnification factor
		zoom = orig_vox_size / vox_size
		
		# When zoom != 1, ie., correct for downward slope of beam		
		dX = get_delta_sx(sx, sx0, zoom)
		dZ = get_delta_sz(tomo_params['energy'], dX)
		
		Z += dZ

		# Positions in units of pixels with correction in case tiles do not have the same number of pixels
		Z = round(Z / vox_size + (Zdim_ref-Zdim)/2.)
		Y = round(Y / vox_size + (Ydim_ref-Ydim)/2.)
		X = round(X / vox_size + (Xdim_ref-Xdim)/2.)
		
		print(X, Y, Z)
	
		path = vol
		line = f"{path} = {X}, {Y}, {Z}"
		line_to_write += line + "\n"
	
	# Write the stitch settings file and copy slurm_config file to local directory
	write_stitch_settings(args.name, args.binning, line_to_write, point_spacing=100, coarse_block_radius=50, coarse_binning=2, cluster_name='SLURM')
	os.rename('./stitch_settings.txt', './' + args.output)
	if not os.path.isfile('./slurm_config.txt'):
		shutil.copyfile('/data/id16a/inhouse1/sware/NRstitcher/pi2-0523/pi2/bin-linux64/release-nocl/slurm_config_id16a.txt', './slurm_config.txt')
	
if __name__ == '__main__':
    main()
	

	
