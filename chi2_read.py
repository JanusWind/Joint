import numpy as np
import matplotlib.pyplot as plt

chi2_file = open( 'chi2_res', 'r' ).readlines()
s = [0]*len( chi2_file )
n = [0]*len( chi2_file )
Rchi2 = [0]*len( chi2_file )

for i in range( len( chi2_file ) ) :
	s[i]     = float(chi2_file[i][0:-1].split(' ')[0])
	n[i]     = float(chi2_file[i][0:-1].split(' ')[1])
	Rchi2[i] = float(chi2_file[i][0:-1].split(' ')[2])

f1, axs1 = plt.subplots( 2, 1, sharex=True, squeeze=True )

axs1[0].axhline( 4.18, c='b', lw='1', label=r'$n_{FC}$')
axs1[0].axhline( 3.99, c='r', lw='1', label=r'$<n_{PL}>$' )
axs1[0].axhline( 3.82, c='r', ls='--', lw='1', label=r'$min(n_{PL})$' )
axs1[0].scatter( s, n, color='k' )
axs1[0].set_ylabel( r'$n$ $(cm^{-3})$', fontsize=18 )
axs1[0].set_xscale( 'log' )
axs1[0].set_ylim( 3.6, 4.4 )
axs1[0].legend(fontsize=18)

#for tick in axs1[0].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[1].scatter( s, Rchi2, color='k' )
axs1[1].set_ylabel( r'${\chi}^2_{FC}/{\chi}^2_{PL}$', fontsize=18 )
axs1[1].set_xscale( 'log' )
axs1[1].set_yscale( 'log' )
axs1[1].axhline( 1, c='g', lw='1' )
axs1[1].set_xlim( 1e-22, 1e-12 )
axs1[1].set_ylim( 1e-5, 1e3 )
axs1[1].set_xlabel( r's $({\chi}^2 \textasciitilde 1/(s{\sigma})^2)$', fontsize=18 )

#for tick in axs1[1].yaxis.get_major_ticks():
#	tick.label.set_fontsize(18)

axs1[0].set_title( r'$T = 1997-01-08/11:59:32$', fontsize=22 )


#plt.xticks( rotation=45, fontsize=16 )

#plt.yticks( fontsize=16 )

plt.subplots_adjust(wspace=0, hspace=0.2)
plt.tight_layout()

#plt.savefig( 'chi-squared_test_plots' +'.eps', bbox_inches='tight', dpi=40 )

