import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
class Project(object):
	# Change EXMPATH to your path to Excimontec before building! 
	EXMPATH = "/Users/benn-m/Documents/Code/Solar_Software/developing_Exm/Excimontec"
	DEFAULT = f"{EXMPATH}/parameters_default.txt"
	EXE = f"{EXMPATH}/Excimontec.exe"
	def __init__(self,project_name,param_dicts,default_name = DEFAULT,overwrite=False):
		"""Updates default parameters with given list
		of paramater dictionaries, with key commented 
		parameter name and single values """
		try:
			f = open(default_name,'r')
			self.default = [line.strip() for line in f.readlines()]
			f.close()
			print('Found Default')
		except FileNotFoundError:
			print(f"Can't find {default_name}")
			raise(FileNotFoundError)
		self.overwrite = overwrite
		self.project_name = project_name
		self.f_names = []
		self.names = [] 
		self.n_sim = len(param_dicts)
		for i,pd in enumerate(param_dicts):
			cf = self.construct_folder(pd,i)

	def construct_folder(self,param_dict,num):
		self.current_params = self.default.copy()
		for key in param_dict.keys():		
			matches = [(ind,line) for ind,line in enumerate(self.default) if key in line]
			if len(matches)==0:
				print(f"Couldn't find parameter {key},check spelling.")
				return 1
			elif len(matches)>1:
				print(f"Found multiple matches for {key}: use a more specific parameter name.")
				return 1
			else:
				ind,line = matches[0]
				val = param_dict[key]
				self.current_params[ind] = self.replace_rc(line,val)
		f_name = self.project_name+'_'+str(num)
		self.f_names.append(f_name)
		try:
			os.mkdir(f_name)
		except FileExistsError:
			if not self.overwrite:
				ow = input(f'Overwrite {f_name}? [Y/n]')
				if ow.lower() != 'y':
					return 2
			shutil.rmtree(f_name)
			os.mkdir(f_name)
		self.write_parameter_file(f'./{f_name}')
		

	def write_parameter_file(self,location):
		f = open(f'{location}/parameters.txt','w+')
		f.write('\n'.join(self.current_params))
		f.close()

	def replace_rc(self,line,val):
		c_ind = line.index('//')
		comment = line[c_ind:]
		return str(val)+' '+comment

	def run(self,exe_loc = EXE,cores=4):
		for f_name in self.f_names:
			# os.chdir(f'./{f_name}')
			os.chdir(f_name)
			print(os.getcwd())
			args = f'mpiexec -n {cores} {exe_loc} ./parameters.txt'
			print(f'{args=}')
			subprocess.run(args,shell=True)
			os.chdir('..')


class Analyzer(object):

	def __init__(self,project_name,param_dicts,f_names):
		self.project_name = project_name
		self.param_dicts = param_dicts
		self.f_names = f_names
		self.load_summaries()
		self.get_results()
		self.n_sims = len(param_dicts)

	def load_summaries(self):
		self.summaries = {} 
		for f_name in self.f_names:
			f = open(f'./{f_name}/analysis_summary.txt',encoding='ISO-8859-1')
			contents = [line.strip()[:-1] for line in f.readlines()]
			self.summaries[f_name] = contents.copy()
			f.close()

	def exciton_diffusion_length(self,summary):
		phrase = 'Exciton diffusion length is '
		matches = [(ind,line) for ind,line in enumerate(summary) if 'Exciton diffusion length is ' in line]
		if len(matches)==0:
			print(f"Couldn't find phrase {phrase},check spelling.")
			return 1
		elif len(matches)>1:
			print(f"Found multiple matches for {phrase}: use a more specific phrase name.")
			return 1
		else:
			ind,line = matches[0]
		is_loc = line.index('is ')
		nm_loc = line.index('nm')
		start = is_loc+3
		end = nm_loc-1
		quantity = line[start:end]
		value = float(quantity[:quantity.index(' ')])
		error = float(quantity[quantity.index(' ')+3:])
		return value,error
		# print(f'{line=}')
		# print(f'{quantity=}')
		# print(f'{value=}')
		# print(f'{error=}')

	def find_value_phrase(self,summary,phrase):
		matches = [(ind,line) for ind,line in enumerate(summary) if phrase in line]
		if len(matches)==0:
			print(f"Couldn't find phrase {phrase},check spelling.")
			return 1
		elif len(matches)>1:
			print(f"Found multiple matches for {phrase}: use a more specific phrase name.")
			return 1
		else:
			ind,line = matches[0]
		is_loc = line.index('is ')
		nm_loc = line.index('nm')
		start = is_loc+3
		end = nm_loc-1
		quantity = line[start:end]
		value = float(quantity[:quantity.index(' ')])
		error = float(quantity[quantity.index(' ')+3:])
		return value,error

	def get_results_csv(self,summary):
		# line = summary.index("CSV formatted results:")
		phrase = "CSV formatted results"
		matches = [(ind,line) for ind,line in enumerate(summary) if phrase in line]
		if len(matches)==0:
			print(f"Couldn't find phrase '{phrase}',check spelling.")
			return 1
		elif len(matches)>1:
			print(f"Found multiple matches for '{phrase}': use a more specific phrase name.")
			return 1
		else:
			ind,line = matches[0]
		names = [name for name in summary[ind+1].split(",")]
		values = np.array([val for val in summary[ind+2].split(",")])
		df = pd.DataFrame(values.reshape(1,15),columns=names)
		# print(f'{values=}')
		# print(f'{names=}')

		return df 

	def get_results(self):
		self.results = []
		for summary in self.summaries.values():
			self.results.append(self.get_results_csv(summary))

	def plot_edl(self,ax):
		values = np.zeros(self.n_sims)
		errors = np.zeros(self.n_sims)
		sigmas = np.zeros(self.n_sims)
		for i,f_name in enumerate(self.f_names):
			summary = self.summaries[f_name]
			values[i],errors[i] = self.exciton_diffusion_length(summary)
			sigmas[i] = self.param_dicts[i]['Energy_stdev_donor']
		ax.errorbar(sigmas,values,yerr=errors,label='Neat \n $R_{0,sh}=3.4*10^{11}$',capsize=10,color='b')
		ax.plot(sigmas,values,"b-o")

		ax.legend()

