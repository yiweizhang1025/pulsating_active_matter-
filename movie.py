from pylab import *
from matplotlib.ticker import  *
from matplotlib.pyplot import  *
from matplotlib import animation
#from matplotlib import cm
from cmocean import cm
from numpy import loadtxt
import time
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
# References
# ----------
# https://matplotlib.org/examples/animation/rain.html
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

# Colors
# ------
# https://nicoguaro.github.io/posts/cyclic_colormaps/
# https://github.com/matplotlib/cmocean
# sudo apt install python-pip
# sudo pip install "cmocean[plots]"


fig_width_pt  = 246.
inches_per_pt = 1./72
fig_width 	  = fig_width_pt*inches_per_pt
fig_size	 	  = [fig_width, fig_width]

params = {
          'text.latex.preamble': [r"\usepackage{amstext}",],
          'axes.linewidth': .5,
          'axes.titley':.2,
          'axes.labelsize': 6.5,
          'text.fontsize':.4,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,	
          'text.usetex': True,
          'figure.figsize': fig_size}
#rcParams.update(params)

left   = 5e-2
bottom = 8e-2
width  = 9e-1
height = 9e-1

# =======================
#  Load data
# =======================

idx       = [1.6]
ep = 1

#param_all = ['movie_adp_rho_%1.1f_omega_1e1_amp_1e0_align_%ge1_lambda_1e-1_trj_1' %(idx[i],ep) for i in range(len(idx))]
param_all = ['movie_adp_rho_%g_omega_1e1_amp_1e0_align_%ge1_cms' %(idx[i],ep) for i in range(len(idx))]
PI = 3.14159265359

N_plot = len(param_all)

# =======================
#  Plot data
# =======================

# Start time counter
start = time.time()

for j in range(N_plot):
	param  = param_all[j]

	# Load parameters
	d_name   = param
	d        = loadtxt(d_name)
	Lx       = d[0]
	Ly       = d[1]
	rho      = d[2]
	T        = d[3]
	omega    = d[6]
	t_total  = d[8]
	t_record = d[10]
	radius   = d[12]
	amp      = d[13]
	lam      = d[15]
    
# Check normalization
# -------------------
#	d = loadtxt(d_name + '_hist')
#	print sum(d[:,1])

	if omega==0.:
		omega = 1
	Np        = int(Lx*Ly*rho)
	nb_record = int(min(1e3,t_total/t_record))

	# Load data
	d_name = param + '_snap'
	count  = 0
	i      = -1
	pos    = zeros((nb_record, Np, 3))
	with open(d_name) as file_name:
		for line in file_name:
			if count%Np==0:
				i  += 1
				if i==nb_record:
					break
			pos[i, count-Np*i, :] = [k for k in line.split()]
			count += 1

	# Plot data
	fig = figure(figsize=[5,5])
	ax  = axes([left, bottom, width, height])
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_axis_off()
	ax.set_title(r"$\rho=%g,\epsilon=%g,\lambda=%g$" %(rho,10*ep,lam),y=.95)
	col = 'dodgerblue'
	alp = .8
	ms  = width*fig.get_figwidth()/(inches_per_pt*Lx)
	L = max(Lx,Ly)
	plot([0, L], [0, L], alpha=0)
	particles = ax.scatter([], [], edgecolor='none', alpha=alp)
	txt       = text((1-2.5e-1)*L/2, -7e-2*L, '', size=1e1)
#	cm = plt.cm.get_cmap('RdYlBu')
	cmap = ["copper", 'RdBu_r', 'Oranges', 'cividis', 'hot', 'plasma']
	# Create a grid with particle size
#	for i in range(int(Lx/radius)):
#		plot([i*radius, i*radius], [0, Ly], 'k', lw=.1)
#	for i in range(int(Ly/radius)):
#		plot([0, Lx], [i*radius, i*radius], 'k', lw=.1)

	# Animate plot
	def animate(i):
			txt.set_text(r'$\omega t/2\pi = %.2g$' %(omega*i*t_record/(2*PI)))
			particles.set_offsets(c_[pos[i,:,0], pos[i,:,1]])
			
#			particles.set_sizes((ones(Np)*radius*ms)**2)
#			particles.set_facecolors(col)
		
			size = radius/(1+lam)*(ones(Np) + lam*sin(pos[i,:,2]))
			particles.set_sizes((size*ms)**2)
			particles.set_facecolors(cm.phase_r(pos[i,:,2]%(2*pi)/(2*pi),1))
#			particles.set_facecolors(cm.hsv(pos[i,:,2]%(2*pi)/(2*pi),1))
			if i == 1: 
				fout = open('lastframe_spiral.dat','w')
				for part in range(Np):
					fout.writelines([str(pos[i,part,0]),'\t',str(pos[i,part,1]),'\t',str(pos[i,part,2]),'\n'])
				fout.close
			print (d_name, ':', i, '/', nb_record)
			return particles, txt

	ani = animation.FuncAnimation(fig, animate, nb_record)#, interval=40, blit=True)
	# Save movie	
	movie_name = 'movie_' + param + '.mp4'
	ani.save(movie_name, fps=15, dpi=3e2)#, bitrate=1e3)

# End script message
duration = time.time() - start
print ('Plot duration:', str('%.3g' %duration), 'seconds')

