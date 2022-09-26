from exm.inter import *
import numpy as np
project_name = "data/exciton_diffusion_parameter_sweep"
sigmas = np.round(np.linspace(0.04,0.10,6),2)
print(sigmas)
folder_names = [f'sigma_{val}' for val in sigmas]
default_parameter_file = './exciton_diffusion_default.txt'
parameter_dicts = [{'Energy_stdev_donor':val,'Energy_stdev_acceptor':val} for val in sigmas]
p = Project(project_name,folder_names,parameter_dicts,default_name=default_parameter_file)
p.run(cores =4)