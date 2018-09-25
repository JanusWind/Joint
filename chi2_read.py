import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams

rcParams['figure.figsize'] = 10, 10

chi2R_file = open( 'chi2R_res', 'r' ).readlines()
s      = [0]*len( chi2R_file )
n      = [0]*len( chi2R_file )
v0_x   = [0]*len( chi2R_file )
v0_y   = [0]*len( chi2R_file )
v0_z   = [0]*len( chi2R_file )
v0     = [0]*len( chi2R_file )
the_v  = [0]*len( chi2R_file )
phi_v  = [0]*len( chi2R_file )
w_per  = [0]*len( chi2R_file )
w_par  = [0]*len( chi2R_file )
w      = [0]*len( chi2R_file )
Rchi2R = [0]*len( chi2R_file )

for i in range( len( chi2R_file ) ) :
	s[i]      = float(chi2R_file[i][0:-1].split(' ')[0])
	n[i]      = float(chi2R_file[i][0:-1].split(' ')[1])
	v0_x[i]   = float(chi2R_file[i][0:-1].split(' ')[2])
	v0_y[i]   = float(chi2R_file[i][0:-1].split(' ')[3])
	v0_z[i]   = float(chi2R_file[i][0:-1].split(' ')[4])
	v0[i]     = float(chi2R_file[i][0:-1].split(' ')[5])
	the_v[i]  = np.rad2deg(np.arccos(  v0_z[i] /  v0[i] ) )
	phi_v[i]  = np.rad2deg(np.arctan2( v0_y[i], v0_x[i] ) )
	w_per[i]  = float(chi2R_file[i][0:-1].split(' ')[6])
	w_par[i]  = float(chi2R_file[i][0:-1].split(' ')[7])
	w[i]      = float(chi2R_file[i][0:-1].split(' ')[8])
	Rchi2R[i] = float(chi2R_file[i][0:-1].split(' ')[9])
#-------FC--------PL
# n     4.18      3.99
# v0_x  -414      -419
# v0_y  33        21
# v0_z  -2        2
# v0    416       419
# the_v 90.275    89.727
# phi_v 175.44    177.13
# w_per --        --
# w_par --        --
# w     25        26

f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )


# n_p

axs1[0].axhline( 4.18, c='b', lw='1', label=r'$n_{FC}$')
axs1[0].axhline( 3.99, c='r', lw='1', label=r'$<n_{PL}>$' )
axs1[0].axhline( 3.82, c='r', ls='--', lw='1', label=r'$min(n_{PL})$' )
axs1[0].scatter( s, n, color='k' )
axs1[0].set_ylabel( r'$n$ $(cm^{-3})$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( 3.6, 4.4 )
axs1[0].legend(fontsize=18)


'''
# v0_x

axs1[0].axhline( -414., c='b', lw='1', label=r'$v_{0xFC}$')
axs1[0].axhline( -419, c='r', lw='1', label=r'$<v_{0xPL}>$' )
axs1[0].scatter( s, v0_x, color='k' )
axs1[0].set_ylabel( r'$v_{0x}$ $(km/s)$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( -412, -421 )
axs1[0].legend(fontsize=18)



# v0_y

axs1[1].axhline( 33., c='b', lw='1', label=r'$v_{0yFC}$')
axs1[1].axhline( 21., c='r', lw='1', label=r'$<v_{0yPL}>$' )
axs1[1].scatter( s, v0_y, color='k' )
axs1[1].set_ylabel( r'$v_{0y}$ $(km/s)$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 20, 34 )
axs1[1].legend(fontsize=18)



# v0_z

axs1[2].axhline( -2., c='b', lw='1', label=r'$v_{0zFC}$')
axs1[2].axhline( 2., c='r', lw='1', label=r'$<v_{0zPL}>$' )
axs1[2].scatter( s, v0_z, color='k' )
axs1[2].set_ylabel( r'$v_{0z}$ $(km/s)$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( -3, 3 )
axs1[2].legend(fontsize=18)


'''
# v0

axs1[1].axhline( 416., c='b', lw='1', label=r'$v_{0FC}$')
axs1[1].axhline( 419., c='r', lw='1', label=r'$<v_{0PL}>$' )
axs1[1].scatter( s, v0, color='k' )
axs1[1].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 415, 420 )
axs1[1].legend(fontsize=18)


'''
# the_v

axs1[1].axhline( 90.275, c='b', lw='1', label=r'$\theta_{vFC}$')
axs1[1].axhline( 89.727, c='r', lw='1', label=r'$<\theta_{vPL}>$' )
axs1[1].scatter( s, the_v, color='k' )
axs1[1].set_ylabel( r'$\theta_{v}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 89.5, 90.5 )
axs1[1].legend(fontsize=18)



# phi_v

axs1[2].axhline( 175.44, c='b', lw='1', label=r'$\phi_{vFC}$')
axs1[2].axhline( 177.13, c='r', lw='1', label=r'$<\phi_{vPL}>$' )
axs1[2].scatter( s, phi_v, color='k' )
axs1[2].set_ylabel( r'$\phi_{v}$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( 175, 178 )
axs1[2].legend(fontsize=18)
'''

'''
# w_per

axs1[1].scatter( s, w_per, color='k' )
axs1[1].set_ylabel( r'$w_{per}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 19, 27 )
axs1[1].legend(fontsize=18)



# w_par

axs1[2].scatter( s, w_par, color='k' )
axs1[2].set_ylabel( r'$w_{par}$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( 19, 27 )
axs1[2].legend(fontsize=18)
'''


# w

axs1[2].axhline( 25, c='b', lw='1', label=r'$w_{FC}$')
axs1[2].axhline( 26, c='r', lw='1', label=r'$<w_{PL}>$' )
axs1[2].scatter( s, w, color='k' )
axs1[2].set_ylabel( r'$w$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( 19, 27 )
axs1[2].legend(fontsize=18)


#for tick in axs1[0].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[3].scatter( s, Rchi2R, color='k' )
axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$', fontsize=18 )
axs1[3].set_xscale( 'log' )
axs1[3].set_yscale( 'log' )
axs1[3].axhline( 1, c='g', lw='1' )
axs1[3].set_xlim( 1e-20, 1e-14 )
axs1[3].set_ylim( 0.1, 10 )
axs1[3].set_xlabel( r's $({\chi}^{2(R)} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

#for tick in axs1[1].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[0].set_title( r'$T = 1997-01-08/11:59:32$', fontsize=22 )


#plt.xticks( rotation=45, fontsize=16 )

#plt.yticks( fontsize=16 )

plt.subplots_adjust(wspace=0, hspace=0.2)
plt.tight_layout()

plt.savefig( 'chi-squared_test_plots_n_v0_w' +'.pdf', bbox_inches='tight', dpi=40 )

plt.show()

