# Creates txt files for importing SPH data into SKIRT
# Be sure to update filePath to the correct location

import pynbody
import numpy as np
from timeit import default_timer as timer
from os.path import expanduser
import argparse

class galaxy:

	parser = argparse.ArgumentParser()
	parser.add_argument("--path")
	args = parser.parse_args()

	# read interval from config.txt (interval = 1 counts every star)
	# store config.txt inputs as dictionary 
	d = {}
	with open(args.path+"/config.txt") as f:
		for line in f:
			(key, val) = line.split()
			d[key] = val

		for name, value in d.items():
			s = name.split('/')

			if s[0] == 'interval':
				interval = int(value)
				print('interval set by config file:', interval)

	# Makes cuts based on position
	xAbsMin = -6.31 * 1e6															
	xAbsMax = -6.24 * 1e6
	yAbsMin = -5.09	* 1e6
	yAbsMax = -5.02	* 1e6
	zAbsMin = 2.98 * 1e6
	zAbsMax = 3.05 * 1e6

	def __init__(self):

		home = expanduser("~")
		filePath = home+"/Research/NIHAO/Data/g8.26e11.01024"

		# Load NIHAO data
		self.data = pynbody.load(filePath)

		self.full_x_pos = np.float32(self.data.star['pos'].in_units('pc')[0:-1:self.interval][:,0])
		self.full_y_pos = np.float32(self.data.star['pos'].in_units('pc')[0:-1:self.interval][:,1])
		self.full_z_pos = np.float32(self.data.star['pos'].in_units('pc')[0:-1:self.interval][:,2])
		self.full_mass = np.float32(self.data.star['massform'].in_units('Msol')[0:-1:self.interval]) 	# in solar masses	
		self.full_metals = np.float32(self.data.star['metals'][0:-1:self.interval])
		self.full_age = np.float32(self.data.star['age'].in_units('yr')[0:-1:self.interval])
		self.full_x_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[0:-1:self.interval][:,0])
		self.full_y_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[0:-1:self.interval][:,1])
		self.full_z_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[0:-1:self.interval][:,2])
		self.full_smooth = np.float32(self.data.star['smooth'].in_units('pc')[0:-1:self.interval])

		self.full_x_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[0:-1:self.interval][:,0])
		self.full_y_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[0:-1:self.interval][:,1])
		self.full_z_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[0:-1:self.interval][:,2])
		self.full_smooth_dust = np.float32(self.data.gas['smooth'].in_units('pc')[0:-1:self.interval])
		self.full_mass_dust = np.float32(self.data.gas['mass'].in_units('Msol')[0:-1:self.interval]) 	# in solar masses	
		self.full_metals_dust = np.float32(self.data.gas['metals'][0:-1:self.interval])
		self.full_temp_dust = np.float32(self.data.gas['temp'][0:-1:self.interval])


		self.full_length_star = len(self.full_x_pos) 						# starting length
		self.starIndex = [] 	
		print(self.full_length_star, 'full star')											# indices to be deleted from data
		
		self.full_length_dust = len(self.full_x_pos_dust) 					# starting length
		self.dustIndex = [] 		
		print(self.full_length_dust, 'full dust')										# indices to be deleted from data
		

	def starCut(self):

		for i in range(self.full_length_star):
			if (self.full_x_pos[i] < self.xAbsMin) or (self.full_x_pos[i] > self.xAbsMax):
				self.starIndex.append(i)
			elif (self.full_y_pos[i] < self.yAbsMin) or (self.full_y_pos[i] > self.yAbsMax):
				self.starIndex.append(i)
			elif (self.full_z_pos[i] < self.zAbsMin) or (self.full_z_pos[i] > self.zAbsMax):
				self.starIndex.append(i)

		self.x_pos = np.float32(np.delete(self.full_x_pos,self.starIndex))
		self.y_pos = np.float32(np.delete(self.full_y_pos,self.starIndex))
		self.z_pos = np.float32(np.delete(self.full_z_pos,self.starIndex))
		self.mass = np.float32(np.delete(self.full_mass,self.starIndex))
		self.metals = np.float32(np.delete(self.full_metals,self.starIndex) )
		self.age = np.float32(np.delete(self.full_age,self.starIndex))
		self.x_vel = np.float32(np.delete(self.full_x_vel,self.starIndex))
		self.y_vel = np.float32(np.delete(self.full_y_vel,self.starIndex))
		self.z_vel = np.float32(np.delete(self.full_z_vel,self.starIndex))
		self.smooth = np.float32(np.delete(self.full_smooth,self.starIndex))


		self.starLength = len(self.x_pos)
		print(self.starLength, 'stars')

	def dustCut(self):

		for i in range(self.full_length_dust):
			if (self.full_x_pos_dust[i] < self.xAbsMin) or (self.full_x_pos_dust[i] > self.xAbsMax):
				self.dustIndex.append(i)
			elif (self.full_y_pos_dust[i] < self.yAbsMin) or (self.full_y_pos_dust[i] > self.yAbsMax):
				self.dustIndex.append(i)
			elif (self.full_z_pos_dust[i] < self.zAbsMin) or (self.full_z_pos_dust[i] > self.zAbsMax):
				self.dustIndex.append(i)

		self.x_pos_dust = np.float32(np.delete(self.full_x_pos_dust,self.dustIndex))
		self.y_pos_dust = np.float32(np.delete(self.full_y_pos_dust,self.dustIndex))
		self.z_pos_dust = np.float32(np.delete(self.full_z_pos_dust,self.dustIndex))
		self.smooth_dust = np.float32(np.delete(self.full_smooth_dust,self.dustIndex))
		self.mass_dust = np.float32(np.delete(self.full_mass_dust,self.dustIndex))
		self.metals_dust = np.float32(np.delete(self.full_metals_dust,self.dustIndex))
		self.temp_dust = np.float32(np.delete(self.full_temp_dust,self.dustIndex))

		self.dustLength = len(self.x_pos_dust)
		print(self.dustLength, 'dust')

	def shift(self):

		xCenter = (self.xAbsMax + self.xAbsMin)/2
		yCenter = (self.yAbsMax + self.yAbsMin)/2
		zCenter = (self.zAbsMax + self.zAbsMin)/2

		self.x_pos = self.x_pos - xCenter
		self.y_pos = self.y_pos - yCenter
		self.z_pos = self.z_pos - zCenter

		self.x_pos_dust = self.x_pos_dust - xCenter
		self.y_pos_dust = self.y_pos_dust - yCenter
		self.z_pos_dust = self.z_pos_dust - zCenter
		

if __name__=='__main__':

	start = timer()

	g = galaxy()

	g.starCut()
	g.dustCut()
	g.shift()

	# compensate for skipping stars
	g.smooth = g.smooth * (g.interval**(1./3.))
	g.mass = g.mass * g.interval

	np.savetxt(g.args.path+'/SKIRT_files/radiation.txt',np.float32(np.c_[g.x_pos, g.y_pos, g.z_pos, g.smooth, g.x_vel, g.y_vel, g.z_vel, g.mass, g.metals, g.age]))
	np.savetxt(g.args.path+'/SKIRT_files/dust.txt',np.float32(np.c_[g.x_pos_dust, g.y_pos_dust, g.z_pos_dust, g.smooth_dust, g.mass_dust, g.metals_dust, g.temp_dust]))

	end = timer()
	print('time: ', end - start)



#a = [1,1,1,1,1]
#b = [2,2,2,2,2]
#c = [3,3,3,3,3]
#np.savetxt('testing.txt',np.c_[a,b,c])