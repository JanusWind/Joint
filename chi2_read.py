import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import os.path

rcParams['figure.figsize'] = 10, 10

#date = '1997-01-09;23:04:27'
date = '1997-01-10;11:59:45'

chi2R_file = open( os.path.join('chi-squared', date,  'chi2R_res_' + date ), 'r' ).readlines()

s        = [0]*len( chi2R_file )
n        = [0]*len( chi2R_file )
n_sig    = [0]*len( chi2R_file )
n_sig_r  = [0]*len( chi2R_file )
gA       = [0]*len( chi2R_file )
gA_sig   = [0]*len( chi2R_file )
gA_sig_r = [0]*len( chi2R_file )
gV       = [0]*len( chi2R_file )
gV_sig   = [0]*len( chi2R_file )
gV_sig_r = [0]*len( chi2R_file )
dthe     = [0]*len( chi2R_file )
dthe_sig = [0]*len( chi2R_file )
dphi     = [0]*len( chi2R_file )
dphi_sig = [0]*len( chi2R_file )
v0_x     = [0]*len( chi2R_file )
v0_y     = [0]*len( chi2R_file )
v0_z     = [0]*len( chi2R_file )
v0       = [0]*len( chi2R_file )
v0_sig   = [0]*len( chi2R_file )
v0_sig_r = [0]*len( chi2R_file )
the_v    = [0]*len( chi2R_file )
phi_v    = [0]*len( chi2R_file )
w_per    = [0]*len( chi2R_file )
w_par    = [0]*len( chi2R_file )
w        = [0]*len( chi2R_file )
Rchi2fc  = [0]*len( chi2R_file )
Rchi2pl  = [0]*len( chi2R_file )
Rchi2R   = [0]*len( chi2R_file )

for i in range( len( chi2R_file ) ) :
	s[i]        = float(chi2R_file[i][0:-1].split(' ')[0])
	n[i]        = float(chi2R_file[i][0:-1].split(' ')[1])
	n_sig[i]    = float(chi2R_file[i][0:-1].split(' ')[2])
	n_sig_r[i]  = n_sig[i] / n[i]
	gA[i]       = float(chi2R_file[i][0:-1].split(' ')[3])
	gA_sig[i]   = float(chi2R_file[i][0:-1].split(' ')[4])
	gA_sig_r[i] = gA_sig[i] / gA[i]
	gV[i]       = float(chi2R_file[i][0:-1].split(' ')[5])
	gV_sig[i]   = float(chi2R_file[i][0:-1].split(' ')[6])
	gV_sig_r[i] = gV_sig[i] / gV[i]
	dthe[i]     = float(chi2R_file[i][0:-1].split(' ')[7])
	dthe_sig[i] = float(chi2R_file[i][0:-1].split(' ')[8])
	dphi[i]     = float(chi2R_file[i][0:-1].split(' ')[9])
	dphi_sig[i] = float(chi2R_file[i][0:-1].split(' ')[10])
	v0_x[i]     = float(chi2R_file[i][0:-1].split(' ')[11])
	v0_y[i]     = float(chi2R_file[i][0:-1].split(' ')[12])
	v0_z[i]     = float(chi2R_file[i][0:-1].split(' ')[13])
	v0[i]       = float(chi2R_file[i][0:-1].split(' ')[14])
	the_v[i]    = np.rad2deg(np.arccos(  v0_z[i] /  v0[i] ) )
	phi_v[i]    = np.rad2deg(np.arctan2( v0_y[i], v0_x[i] ) )
	w_per[i]    = float(chi2R_file[i][0:-1].split(' ')[15])
	w_par[i]    = float(chi2R_file[i][0:-1].split(' ')[16])
	w[i]        = float(chi2R_file[i][0:-1].split(' ')[17])
	Rchi2fc[i]  = float(chi2R_file[i][0:-1].split(' ')[18])
	Rchi2pl[i]  = float(chi2R_file[i][0:-1].split(' ')[19])
	Rchi2R[i]   = Rchi2fc[i]/Rchi2pl[i]

#-------------------------#
#---1997-01-09/23:04:27---#
#-------FC--------PL------#
# n     8.35      8.4
# v0_x  -373      -379
# v0_y  30        18
# v0_z  7         11
# v0    374       380
# the_v 88.928    88.339
# phi_v 175.40    177.28
# w_per --        --
# w_par --        --
# w     18        21

best_gA = min( gA_sig )
best_i = np.where( np.array( gA_sig ) == best_gA )[0][0]
best_s = s[best_i]


print 'data loaded'

def plot_all() :

	plot_n_gA_sig()

	plot_v0_gV_sig()

	plot_n_v0_w()

	plot_v0xyz()

	plot_v0_the_phi()

	plot_the_dthe_phi_dphi()

	plot_w_per_par()

	plot_chi2()

	plot_gA_gV_dthe_dphi()



def plot_n_gA_sig() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# n_p_c

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( n[best_i], c='b', lw='1', label=r'$n(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(n[best_i], 4) ) )
	axs1[0].scatter( s, n, color='k' )
	axs1[0].errorbar( s, n, yerr=n_sig, color='k', ls='', elinewidth=2 )
	axs1[0].set_ylabel( r'$n_{pc}$ $(cm^{-3})$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# n_p_c_sig

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( n_sig[best_i], c='b', lw='1', label=r'$\sigma_n(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(n_sig[best_i], 4) ) )
	axs1[1].scatter( s, n_sig, color='k' )
	axs1[1].set_ylabel( r'$\sigma_n$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_yscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	# gA

	axs1[2].axvline( best_s, c='r', lw='1' )
	axs1[2].axhline( gA[best_i], c='b', lw='1', label=r'$g_A(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gA[best_i], 4) ) )
	axs1[2].scatter( s, gA, color='k' )
	axs1[2].errorbar( s, gA, yerr=gA_sig, color='k', ls='', elinewidth=2 )
	axs1[2].set_ylabel( r'$g_A$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].legend(fontsize=18)


	#gA_sig

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( min(gA_sig), color='r', lw='1', label=r'$min(\sigma_{g_A})=$'+str(round(min(gA_sig), 4) ) )
	axs1[3].scatter( s, gA_sig, color='k' )
	axs1[3].set_ylabel( r'$\sigma_{g_A}$', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_yscale( 'log' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].legend(fontsize=18)

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'n_gA_sig'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()




def plot_v0_gV_sig() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# v0

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( v0[best_i], c='b', lw='1', label=r'$v_0(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0[best_i], 4) ) )
	axs1[0].scatter( s, v0, color='k' )
	axs1[0].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# gV

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( gV[best_i], c='b', lw='1', label=r'$g_V(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV[best_i], 4) ) )
	axs1[1].scatter( s, gV, color='k' )
	axs1[1].errorbar( s, gV, yerr=gV_sig, color='k', ls='', elinewidth=2 )
	axs1[1].set_ylabel( r'$g_V$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	#gV_sig

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( gV_sig[best_i], c='b', lw='1', label=r'$\sigma_{g_V}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV_sig[best_i], 4) ) )
	axs1[2].axhline( min(gV_sig), color='r', lw='1', ls='--', label=r'$min(\sigma_{g_V})=$'+str(round(min(gV_sig), 6) ) )
	axs1[2].scatter( s, gV_sig, color='k' )
	axs1[2].set_ylabel( r'$\sigma_{g_V}$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_yscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].legend(fontsize=18)

	# Reduced chi-squared Ratio

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( Rchi2R[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2R[best_i], 4) ) )
	axs1[3].scatter( s, Rchi2R, color='k' )
	axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_yscale( 'log' )
	axs1[3].axhline( 1, c='g', lw='1' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].set_ylim( 1e14, 1e19 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'v0_gV_sig'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()






def plot_n_v0_w() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# n_p_c

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( n[best_i], c='b', lw='1', label=r'$n(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(n[best_i], 4) ) )
	axs1[0].scatter( s, n, color='k' )
	axs1[0].errorbar( s, n, yerr=n_sig, color='k', ls='', elinewidth=2 )
	axs1[0].set_ylabel( r'$n_{pc}$ $(cm^{-3})$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# v0

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( v0[best_i], c='b', lw='1', label=r'$v_0(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0[best_i], 4) ) )
	axs1[1].scatter( s, v0, color='k' )
	axs1[1].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	# w

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( w[best_i], c='b', lw='1', label=r'$w(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(w[best_i], 4) ) )
	axs1[2].scatter( s, w, color='k' )
	axs1[2].set_ylabel( r'$w$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].legend(fontsize=18)

	# Reduced chi-squared Ratio

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( Rchi2R[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2R[best_i], 4) ) )
	axs1[3].scatter( s, Rchi2R, color='k' )
	axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_yscale( 'log' )
	axs1[3].axhline( 1, c='g', lw='1' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].set_ylim( 1e14, 1e19 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'n_v0_w'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()




def plot_v0xyz() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# v0

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( v0[best_i], c='b', lw='1', label=r'$v_0(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0[best_i], 4) ) )
	axs1[0].scatter( s, v0, color='k' )
	axs1[0].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# v0_x

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( v0_x[best_i], c='b', lw='1', label=r'$v_{0x}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0_x[best_i], 4) ) )
	axs1[1].scatter( s, v0_x, color='k' )
	axs1[1].set_ylabel( r'$v_{0x}$ $(km/s)$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	# v0_y

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( v0_y[best_i], c='b', lw='1', label=r'$v_{0y}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0_y[best_i], 4) ) )
	axs1[2].scatter( s, v0_y, color='k' )
	axs1[2].set_ylabel( r'$v_{0y}$ $(km/s)$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].legend(fontsize=18)

	# v0_z

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( v0_z[best_i], c='b', lw='1', label=r'$v_{0z}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0_z[best_i], 4) ) )
	axs1[3].scatter( s, v0_z, color='k' )
	axs1[3].set_ylabel( r'$v_{0z}$ $(km/s)$', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'v0xyz'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()




def plot_v0_the_phi() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# v0

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( v0[best_i], c='b', lw='1', label=r'$v_0(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(v0[best_i], 4) ) )
	axs1[0].scatter( s, v0, color='k' )
	axs1[0].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# the_v

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( the_v[best_i], c='b', lw='1', label=r'$\theta_v(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(the_v[best_i], 4) ) )
	axs1[1].scatter( s, the_v, color='k' )
	axs1[1].set_ylabel( r'$\theta_{v}$ $(^\circ)$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	# phi_v

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( phi_v[best_i], c='b', lw='1', label=r'$\phi_v(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(phi_v[best_i], 4) ) )
	axs1[2].scatter( s, phi_v, color='k' )
	axs1[2].set_ylabel( r'$\phi_{v}$ $(^\circ)$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].legend(fontsize=18)

	# Reduced chi-squared Ratio

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( Rchi2R[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2R[best_i], 4) ) )
	axs1[3].scatter( s, Rchi2R, color='k' )
	axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_yscale( 'log' )
	axs1[3].axhline( 1, c='g', lw='1' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].set_ylim( 1e14, 1e19 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'v0_the_phi'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()



def plot_the_dthe_phi_dphi() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# the_v

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( the_v[best_i], c='b', lw='1', label=r'$\theta_v(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(the_v[best_i], 4) ) )
	axs1[0].scatter( s, the_v, color='k' )
	axs1[0].set_ylabel( r'$\theta_{v}$ $(^\circ)$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	#axs1[0].set_ylim( 88.3, 89. )
	axs1[0].legend(fontsize=18)

	# dthe

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( dthe[best_i], c='b', lw='1', label=r'$\Delta\theta(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dthe[best_i], 4) ) )
	axs1[1].scatter( s, dthe, color='k' )
	axs1[1].set_ylabel( r'$\Delta\theta$ $(^\circ)$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].errorbar( s, dthe, yerr=dthe_sig, color='k', ls='', elinewidth=2 )
	axs1[1].legend(fontsize=18)

	# phi_v

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( phi_v[best_i], c='b', lw='1', label=r'$\phi_v(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(phi_v[best_i], 4) ) )
	axs1[2].scatter( s, phi_v, color='k' )
	axs1[2].set_ylabel( r'$\phi_{v}$ $(^\circ)$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	#axs1[2].set_ylim( 175, 178 )
	axs1[2].legend(fontsize=18)

	# dphi

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( dphi[best_i], c='b', lw='1', label=r'$\Delta\phi(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dphi[best_i], 4) ) )
	axs1[3].scatter( s, dphi, color='k' )
	axs1[3].set_ylabel( r'$\Delta\phi$ $(^\circ)$', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].errorbar( s, dphi, yerr=dphi_sig, color='k', ls='', elinewidth=2 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'the_dthe_phi_dphi'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()





def plot_w_per_par() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# w

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( w[best_i], c='b', lw='1', label=r'$w(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(w[best_i], 4) ) )
	axs1[0].scatter( s, w, color='k' )
	axs1[0].set_ylabel( r'$w$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	#axs1[0].set_ylim( 13, 22 )
	axs1[0].legend(fontsize=18)

	# w_per

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( w_per[best_i], c='b', lw='1', label=r'$w_{\perp}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(w_per[best_i], 4) ) )
	axs1[1].scatter( s, w_per, color='k' )
	axs1[1].set_ylabel( r'$w_{\perp}$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	#axs1[1].set_ylim( 12, 22 )
	axs1[1].legend(fontsize=18)

	# w_par

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( w_par[best_i], c='b', lw='1', label=r'$w_{\parallel}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(w_par[best_i], 4) ) )
	axs1[2].scatter( s, w_par, color='k' )
	axs1[2].set_ylabel( r'$w_{\parallel}$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	#axs1[2].set_ylim( 14, 22 )
	axs1[2].legend(fontsize=18)

	# Reduced chi-squared Ratio

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( Rchi2R[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2R[best_i], 4) ) )
	axs1[3].scatter( s, Rchi2R, color='k' )
	axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_yscale( 'log' )
	axs1[3].axhline( 1, c='g', lw='1' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].set_ylim( 1e14, 1e19 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'w_per_par'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()




def plot_chi2() :

	f1, axs1 = plt.subplots( 3, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# Reduced chi-squared FC

	axs1[0].axvline( best_s, c='r', lw='1' ) 
	axs1[0].axhline( Rchi2fc[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2fc[best_i], 4) ) )
	axs1[0].scatter( s, Rchi2fc, color='k' )
	axs1[0].set_ylabel( r'${\chi}^{2(R)}_{FC}$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_yscale( 'log' )
	#axs1[0].axhline( 1, c='g', lw='1' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# Reduced chi-squared PL

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( Rchi2pl[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(Rchi2pl[best_i] ) )
	axs1[1].scatter( s, Rchi2pl, color='k' )
	axs1[1].set_ylabel( r'${\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_yscale( 'log' )
	axs1[1].axhline( 1, c='g', lw='1' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].set_ylim( 1e-13, 1e-12 )
	axs1[1].legend(fontsize=18)

	# Reduced chi-squared Ratio


	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( Rchi2R[best_i], c='b', lw='1', label=r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(Rchi2R[best_i], 4) ) )
	axs1[2].scatter( s, Rchi2R, color='k' )
	axs1[2].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_yscale( 'log' )
	axs1[2].axhline( 1, c='g', lw='1' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].set_ylim( 1e14, 1e19 )
	axs1[2].legend(fontsize=18)

	axs1[2].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'chi2'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()


def plot_gA_gV_dthe_dphi() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# gA

	axs1[0].axvline( best_s, c='r', lw='1' )
	axs1[0].axhline( gA[best_i], c='b', lw='1', label=r'$g_A(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gA[best_i], 4) ) )
	axs1[0].scatter( s, gA, color='k' )
	axs1[0].errorbar( s, gA, yerr=gA_sig, color='k', ls='', elinewidth=2 )
	axs1[0].set_ylabel( r'$g_A$', fontsize=18 )
	axs1[0].set_xscale( 'log' )
	axs1[0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0].legend(fontsize=18)

	# gV

	axs1[1].axvline( best_s, c='r', lw='1' ) 
	axs1[1].axhline( gV[best_i], c='b', lw='1', label=r'$g_V(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV[best_i], 4) ) )
	axs1[1].scatter( s, gV, color='k' )
	axs1[1].errorbar( s, gV, yerr=gV_sig, color='k', ls='', elinewidth=2 )
	axs1[1].set_ylabel( r'$g_V$', fontsize=18 )
	axs1[1].set_xscale( 'log' )
	axs1[1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1].legend(fontsize=18)

	# dthe

	axs1[2].axvline( best_s, c='r', lw='1' ) 
	axs1[2].axhline( dthe[best_i], c='b', lw='1', label=r'$\Delta\theta(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dthe[best_i], 4) ) )
	axs1[2].scatter( s, dthe, color='k' )
	axs1[2].set_ylabel( r'$\Delta\theta$ $(^\circ)$', fontsize=18 )
	axs1[2].set_xscale( 'log' )
	axs1[2].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2].errorbar( s, dthe, yerr=dthe_sig, color='k', ls='', elinewidth=2 )
	axs1[2].legend(fontsize=18)

	# dphi

	axs1[3].axvline( best_s, c='r', lw='1' ) 
	axs1[3].axhline( dphi[best_i], c='b', lw='1', label=r'$\Delta\phi(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dphi[best_i], 4) ) )
	axs1[3].scatter( s, dphi, color='k' )
	axs1[3].set_ylabel( r'$\Delta\phi$ $(^\circ)$', fontsize=18 )
	axs1[3].set_xscale( 'log' )
	axs1[3].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3].errorbar( s, dphi, yerr=dphi_sig, color='k', ls='', elinewidth=2 )
	axs1[3].legend(fontsize=18)

	axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'gA_gV_dthe_dphi'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()

def plot_gA_gV_dthe_dphi_sig() :

	f1, axs1 = plt.subplots( 4, 2, sharex=True, squeeze=True )

	axs1[0][0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# gA

	axs1[0][0].axvline( best_s, c='r', lw='1' )
	axs1[0][0].axhline( gA[best_i], c='b', lw='1', label=r'$g_A(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gA[best_i], 4) ) )
	axs1[0][0].scatter( s, gA, color='k' )
	axs1[0][0].errorbar( s, gA, yerr=gA_sig, color='k', ls='', elinewidth=2 )
	axs1[0][0].set_ylabel( r'$g_A$', fontsize=18 )
	axs1[0][0].set_xscale( 'log' )
	axs1[0][0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0][0].legend(fontsize=14)

	#gA_sig

	axs1[1][0].axvline( best_s, c='r', lw='1' ) 
	axs1[1][0].axhline( min(gA_sig), color='b', lw='1', label=r'$min(\sigma_{g_A})=$'+str(round(min(gA_sig), 4) ) )
	axs1[1][0].scatter( s, gA_sig, color='k' )
	axs1[1][0].set_ylabel( r'$\sigma_{g_A}$', fontsize=18 )
	axs1[1][0].set_xscale( 'log' )
	axs1[1][0].set_yscale( 'log' )
	axs1[1][0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1][0].legend(fontsize=14)

	# gV

	axs1[2][0].axvline( best_s, c='r', lw='1' ) 
	axs1[2][0].axhline( gV[best_i], c='b', lw='1', label=r'$g_V(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV[best_i], 4) ) )
	axs1[2][0].scatter( s, gV, color='k' )
	axs1[2][0].errorbar( s, gV, yerr=gV_sig, color='k', ls='', elinewidth=2 )
	axs1[2][0].set_ylabel( r'$g_V$', fontsize=18 )
	axs1[2][0].set_xscale( 'log' )
	axs1[2][0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2][0].legend(fontsize=14)

	#gV_sig

	axs1[3][0].axvline( best_s, c='r', lw='1' ) 
	axs1[3][0].axhline( gV_sig[best_i], c='b', lw='1', label=r'$\sigma_{g_V}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV_sig[best_i], 4) ) )
	axs1[3][0].axhline( min(gV_sig), color='r', lw='1', ls='--', label=r'$min(\sigma_{g_V})=$'+str(round(min(gV_sig), 6) ) )
	axs1[3][0].scatter( s, gV_sig, color='k' )
	axs1[3][0].set_ylabel( r'$\sigma_{g_V}$', fontsize=18 )
	axs1[3][0].set_xscale( 'log' )
	axs1[3][0].set_yscale( 'log' )
	axs1[3][0].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3][0].legend(fontsize=14)

	axs1[3][0].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	axs1[0][1].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )

	# dthe

	axs1[0][1].axvline( best_s, c='r', lw='1' ) 
	axs1[0][1].axhline( dthe[best_i], c='b', lw='1', label=r'$\Delta\theta(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dthe[best_i], 4) ) )
	axs1[0][1].scatter( s, dthe, color='k' )
	axs1[0][1].set_ylabel( r'$\Delta\theta$ $(^\circ)$', fontsize=18 )
	axs1[0][1].set_xscale( 'log' )
	axs1[0][1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[0][1].errorbar( s, dthe, yerr=dthe_sig, color='k', ls='', elinewidth=2 )
	axs1[0][1].legend(fontsize=14)

	# dthe_sig

	axs1[1][1].axvline( best_s, c='r', lw='1' ) 
	axs1[1][1].axhline( dthe_sig[best_i], c='b', lw='1', label=r'$\sigma_{\Delta\theta}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dthe_sig[best_i], 4) ) )
	axs1[1][1].axhline( min(dthe_sig), color='r', lw='1', ls='--', label=r'$min(\sigma_{\Delta\theta})=$'+str(round(min(dthe_sig), 4) ) )
	axs1[1][1].scatter( s, dthe_sig, color='k' )
	axs1[1][1].set_ylabel( r'$\sigma_{\Delta\theta}$ $(^\circ)$', fontsize=18 )
	axs1[1][1].set_xscale( 'log' )
	axs1[1][1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[1][1].legend(fontsize=14)

	# dphi

	axs1[2][1].axvline( best_s, c='r', lw='1' ) 
	axs1[2][1].axhline( dphi[best_i], c='b', lw='1', label=r'$\Delta\phi(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dphi[best_i], 4) ) )
	axs1[2][1].scatter( s, dphi, color='k' )
	axs1[2][1].set_ylabel( r'$\Delta\phi$ $(^\circ)$', fontsize=18 )
	axs1[2][1].set_xscale( 'log' )
	axs1[2][1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[2][1].errorbar( s, dphi, yerr=dphi_sig, color='k', ls='', elinewidth=2 )
	axs1[2][1].legend(fontsize=14)

	# dphi_sig

	axs1[3][1].axvline( best_s, c='r', lw='1' ) 
	axs1[3][1].axhline( dphi_sig[best_i], c='b', lw='1', label=r'$\sigma_{\Delta\phi}(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dphi_sig[best_i], 4) ) )
	axs1[3][1].axhline( min(dphi_sig), color='r', lw='1', ls='--', label=r'$min(\sigma_{\Delta\phi})=$'+str(round(min(dphi_sig), 4) ) )
	axs1[3][1].scatter( s, dphi_sig, color='k' )
	axs1[3][1].set_ylabel( r'$\sigma_{\Delta\phi}$ $(^\circ)$', fontsize=18 )
	axs1[3][1].set_xscale( 'log' )
	axs1[3][1].set_xlim( 5.62e-19, 1.778e-15 )
	axs1[3][1].legend(fontsize=14)

	axs1[3][1].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'gA_gV_dthe_dphi_sig'

	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()

'''
# v0_sig

axs1[1].axhline( min(v0_sig), color='r', lw='1', label=r'$min(\sigma_{v0})=${}'.format(round(min(v0_sig), 6) ) )
axs1[1].scatter( s, v0_sig, color='k' )
axs1[1].set_ylabel( r'$\sigma_{v0}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_yscale( 'log' )
axs1[1].set_xlim( 9e-21, 2e-14 )
#axs1[0].set_ylim( 7.3, 8.5 )
axs1[1].legend(fontsize=18)
'''
