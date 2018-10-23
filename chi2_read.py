import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import os.path

rcParams['figure.figsize'] = 10, 10

date = '1997-01-09;23:04:27'

chi2R_file = open( os.path.join('chi-squared', date,  'chi2R_res_' + date ), 'r' ).readlines()

s        = [0]*len( chi2R_file )
n        = [0]*len( chi2R_file )
n_sig    = [0]*len( chi2R_file )
n_sig_r  = [0]*len( chi2R_file )
gn       = [0]*len( chi2R_file )
gn_sig   = [0]*len( chi2R_file )
gn_sig_r = [0]*len( chi2R_file )
gV       = [0]*len( chi2R_file )
gV_sig   = [0]*len( chi2R_file )
gV_sig_r = [0]*len( chi2R_file )
v0_x     = [0]*len( chi2R_file )
v0_y     = [0]*len( chi2R_file )
v0_z     = [0]*len( chi2R_file )
v0       = [0]*len( chi2R_file )
the_v    = [0]*len( chi2R_file )
phi_v    = [0]*len( chi2R_file )
w_per    = [0]*len( chi2R_file )
w_par    = [0]*len( chi2R_file )
w        = [0]*len( chi2R_file )
Rchi2fc  = [0]*len( chi2R_file )
Rchi2pl  = [0]*len( chi2R_file )
Rchi2R   = [0]*len( chi2R_file )

for i in range( len( chi2R_file ) ) :
	s[i]       = float(chi2R_file[i][0:-1].split(' ')[0])
	n[i]       = float(chi2R_file[i][0:-1].split(' ')[1])
	n_sig[i]   = float(chi2R_file[i][0:-1].split(' ')[2])
	n_sig_r[i] = n_sig[i] / n[i]
	gn[i]       = float(chi2R_file[i][0:-1].split(' ')[3])
	gn_sig[i]   = float(chi2R_file[i][0:-1].split(' ')[4])
	gn_sig_r[i] = gn_sig[i] / gn[i]
	gV[i]       = float(chi2R_file[i][0:-1].split(' ')[5])
	gV_sig[i]   = float(chi2R_file[i][0:-1].split(' ')[6])
	gV_sig_r[i] = gV_sig[i] / gV[i]
	v0_x[i]    = float(chi2R_file[i][0:-1].split(' ')[7])
	v0_y[i]    = float(chi2R_file[i][0:-1].split(' ')[8])
	v0_z[i]    = float(chi2R_file[i][0:-1].split(' ')[9])
	v0[i]      = float(chi2R_file[i][0:-1].split(' ')[10])
	the_v[i]   = np.rad2deg(np.arccos(  v0_z[i] /  v0[i] ) )
	phi_v[i]   = np.rad2deg(np.arctan2( v0_y[i], v0_x[i] ) )
	w_per[i]   = float(chi2R_file[i][0:-1].split(' ')[11])
	w_par[i]   = float(chi2R_file[i][0:-1].split(' ')[12])
	w[i]       = float(chi2R_file[i][0:-1].split(' ')[13])
	Rchi2fc[i] = float(chi2R_file[i][0:-1].split(' ')[14])
	Rchi2pl[i] = float(chi2R_file[i][0:-1].split(' ')[15])
	Rchi2R[i]  = Rchi2fc[i]/Rchi2pl[i]

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

f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )
'''
# n_p_c

axs1[0].axhline( 8.35, c='b', lw='1', label=r'$n_{FC}$' )
axs1[0].axhline( 8.4, c='r', lw='1', label=r'$<n_{PL}>$' )
axs1[0].axhline( 8.05, c='r', ls='--', lw='1', label=r'$min(n_{PL})$' )
axs1[0].scatter( s, n, color='k' )
axs1[0].errorbar( s, n, yerr=n_sig, color='k', ls='', elinewidth=2 )
axs1[0].set_ylabel( r'$n_{pc}$ $(cm^{-3})$, $\sigma_n$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_xlim( 9e-21, 2e-14 )
axs1[0].set_ylim( 8, 10 )
axs1[0].legend(fontsize=18)
'''
'''
# n_p_c_err_r

axs1[1].axhline( min(n_sig_r), color='r', lw='1', label=r'$min(\sigma_n/n)=${}'.format(round(min(n_sig_r), 6) ) )
axs1[1].scatter( s, n_sig_r, color='k' )
axs1[1].set_ylabel( r'$\sigma_n/n$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_yscale( 'log' )
axs1[1].set_xlim( 9e-21, 2e-14 )
#axs1[0].set_ylim( 7.3, 8.5 )
axs1[1].legend(fontsize=18)

'''
'''
# gn

axs1[1].axhline( max(gn), c='b', lw='1', label=r'$max(g_n)$')
axs1[1].axhline( min(gn), c='r', lw='1', label=r'$min(g_n)$' )
axs1[1].scatter( s, gn, color='k' )
axs1[1].errorbar( s, gn, yerr=gn_sig, color='k', ls='', elinewidth=2 )
axs1[1].set_ylabel( r'$g_n$, $\sigma_{g_n}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_xlim( 9e-21, 2e-14 )
axs1[1].set_ylim( 0.5, 1.1 )
axs1[1].legend(fontsize=18)
'''
'''

#gn_err_r

axs1[3].axhline( min(gn_sig_r), color='r', lw='1', label=r'$min(\sigma_{g_n}/{g_n})=$'+str(round(min(gn_sig_r), 6) ) )
axs1[3].scatter( s, gn_sig_r, color='k' )
axs1[3].set_ylabel( r'$\sigma_{g_n}/g_n$', fontsize=18 )
axs1[3].set_xscale( 'log' )
axs1[3].set_yscale( 'log' )
axs1[3].set_xlim( 9e-21, 2e-14 )
#axs1[1].set_ylim( 0.5, 1.1 )
axs1[3].legend(fontsize=18)
'''


'''
# v0_x

axs1[0].axhline( -373., c='b', lw='1', label=r'$v_{0xFC}$')
axs1[0].axhline( -379, c='r', lw='1', label=r'$<v_{0xPL}>$' )
axs1[0].scatter( s, v0_x, color='k' )
axs1[0].set_ylabel( r'$v_{0x}$ $(km/s)$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( -380, -372 )
axs1[0].legend(fontsize=18)



# v0_y

axs1[1].axhline( 30., c='b', lw='1', label=r'$v_{0yFC}$')
axs1[1].axhline( 18., c='r', lw='1', label=r'$<v_{0yPL}>$' )
axs1[1].scatter( s, v0_y, color='k' )
axs1[1].set_ylabel( r'$v_{0y}$ $(km/s)$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 17, 31 )
axs1[1].legend(fontsize=18)



# v0_z

axs1[2].axhline( 7., c='b', lw='1', label=r'$v_{0zFC}$')
axs1[2].axhline( 11., c='r', lw='1', label=r'$<v_{0zPL}>$' )
axs1[2].scatter( s, v0_z, color='k' )
axs1[2].set_ylabel( r'$v_{0z}$ $(km/s)$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( 6, 12 )
axs1[2].legend(fontsize=18)
'''


# v0

axs1[0].axhline( 374., c='b', lw='1', label=r'$v_{0FC}$')
axs1[0].axhline( 380., c='r', lw='1', label=r'$<v_{0PL}>$' )
axs1[0].scatter( s, v0, color='k' )
axs1[0].set_ylabel( r'$v_{0}$ $(km/s)$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( 373, 381 )
axs1[0].legend(fontsize=18)


'''
# v0_err_r

axs1[1].axhline( min(v0_sig_r), color='r', lw='1', label=r'$min(\sigma_{v0}/v0)=${}'.format(round(min(v0_sig_r), 6) ) )
axs1[1].scatter( s, v0_sig_r, color='k' )
axs1[1].set_ylabel( r'$\sigma_{v0}/v0$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_yscale( 'log' )
axs1[1].set_xlim( 9e-21, 2e-14 )
#axs1[0].set_ylim( 7.3, 8.5 )
axs1[1].legend(fontsize=18)
'''

# gV

axs1[2].axhline( max(gV), c='b', lw='1', label=r'$max(g_V)$')
axs1[2].axhline( min(gV), c='r', lw='1', label=r'$min(g_V)$' )
axs1[2].scatter( s, gV, color='k' )
axs1[2].errorbar( s, gV, yerr=gV_sig, color='k', ls='', elinewidth=2 )
axs1[2].set_ylabel( r'$g_V$, $\sigma_{g_V}$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_xlim( 9e-21, 2e-14 )
#axs1[2].set_ylim( 0.5, 1.1 )
axs1[2].legend(fontsize=18)


#gV_err_r

axs1[3].axhline( min(gV_sig_r), color='r', lw='1', label=r'$min(\sigma_{g_V}/{g_V})=$'+str(round(min(gV_sig_r), 6) ) )
axs1[3].scatter( s, gV_sig_r, color='k' )
axs1[3].set_ylabel( r'$\sigma_{g_V}/g_V$', fontsize=18 )
axs1[3].set_xscale( 'log' )
axs1[3].set_yscale( 'log' )
axs1[3].set_xlim( 9e-21, 2e-14 )
#axs1[1].set_ylim( 0.5, 1.1 )
axs1[3].legend(fontsize=18)





'''
# the_v

axs1[1].axhline( 88.928, c='b', lw='1', label=r'$\theta_{vFC}$')
axs1[1].axhline( 88.339, c='r', lw='1', label=r'$<\theta_{vPL}>$' )
axs1[1].scatter( s, the_v, color='k' )
axs1[1].set_ylabel( r'$\theta_{v}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_ylim( 88.3, 89. )
axs1[1].legend(fontsize=18)



# phi_v

axs1[2].axhline( 175.40, c='b', lw='1', label=r'$\phi_{vFC}$')
axs1[2].axhline( 177.28, c='r', lw='1', label=r'$<\phi_{vPL}>$' )
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
axs1[1].set_ylim( 12, 22 )
axs1[1].legend(fontsize=18)



# w_par

axs1[2].scatter( s, w_par, color='k' )
axs1[2].set_ylabel( r'$w_{par}$', fontsize=18 )
axs1[2].set_xscale( 'log' )
axs1[2].set_ylim( 14, 22 )
axs1[2].legend(fontsize=18)



# w

axs1[0].axhline( 18, c='b', lw='1', label=r'$w_{FC}$')
axs1[0].axhline( 21, c='r', lw='1', label=r'$<w_{PL}>$' )
axs1[0].scatter( s, w, color='k' )
axs1[0].set_ylabel( r'$w$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( 13, 22 )
axs1[0].legend(fontsize=18)
'''


'''
#for tick in axs1[0].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[0].scatter( s, Rchi2fc, color='k' )
axs1[0].set_ylabel( r'${\chi}^{2(R)}_{FC}$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_yscale( 'log' )
#axs1[0].axhline( 1, c='g', lw='1' )
axs1[0].set_xlim( 1e-20, 1e-14 )
#axs1[0].set_ylim( 0.01, 100 )

axs1[1].scatter( s, Rchi2pl, color='k' )
axs1[1].set_ylabel( r'${\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_yscale( 'log' )
axs1[1].axhline( 1, c='g', lw='1' )
axs1[1].set_xlim( 1e-20, 1e-14 )
axs1[1].set_ylim( 2e-13, 4e-13 )
#axs1[1].set_ylim( 0.01, 100 )
'''
'''
axs1[3].scatter( s, Rchi2R, color='k' )
axs1[3].set_ylabel( r'${\chi}^{2(R)}_{FC}/{\chi}^{2(R)}_{PL}$ (no s)', fontsize=18 )
axs1[3].set_xscale( 'log' )
axs1[3].set_yscale( 'log' )
axs1[3].axhline( 1, c='g', lw='1' )
axs1[3].set_xlim( 1e-20, 1e-14 )
axs1[3].set_ylim( 1e14, 1e17 )
#axs1[2].set_ylim( 0.01, 100 )
'''
axs1[3].set_xlabel( r's $({\chi}^{2(R)}_{PL} \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

#for tick in axs1[1].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[0].set_title( r'$T = $'+date.replace(';','/'), fontsize=22 )


#plt.xticks( rotation=45, fontsize=16 )

#plt.yticks( fontsize=16 )

plt.subplots_adjust(wspace=0, hspace=0.2)
plt.tight_layout()

plot_vars = 'v0_gV_sig_r'

plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
plt.savefig( os.path.join( 'chi-squared', date, plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

plt.show()

