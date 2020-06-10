import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import os.path
from datetime import datetime, timedelta

rcParams['figure.figsize'] = 10, 10

dates = []

dates += ['1997-01-09;03:01:01']
dates += ['1997-01-09;03:10:42']
dates += ['1997-01-09;03:21:37']
dates += ['1997-01-09;03:31:14']
dates += ['1997-01-09;03:40:52']
dates += ['1997-01-09;03:50:27']
dates += ['1997-01-09;04:00:04']
dates += ['1997-01-09;04:09:42']
dates += ['1997-01-09;04:20:36']
dates += ['1997-01-09;04:30:17']
dates += ['1997-01-09;04:39:49']
dates += ['1997-01-09;04:49:26']
dates += ['1997-01-09;05:00:27']
dates += ['1997-01-09;05:09:58']
dates += ['1997-01-09;05:19:39']
dates += ['1997-01-09;05:30:37']
dates += ['1997-01-09;05:40:11']
dates += ['1997-01-09;05:49:52']
dates += ['1997-01-09;06:03:29']
dates += ['1997-01-09;06:10:24']
dates += ['1997-01-09;06:20:02']
dates += ['1997-01-09;06:29:33']
dates += ['1997-01-09;06:40:37']
dates += ['1997-01-09;06:48:49']

prms_file = open( os.path.join('chi-squared', 'intercal_params_' + dates[0] + '_' + dates[-1] ), 'r' ).readlines()

date_str = [0]*len( prms_file )
date     = [0]*len( prms_file )
s        = [0]*len( prms_file )
gA       = [0]*len( prms_file )
gA_sig   = [0]*len( prms_file )
gV       = [0]*len( prms_file )
gV_sig   = [0]*len( prms_file )
dthe     = [0]*len( prms_file )
dthe_sig = [0]*len( prms_file )
dphi     = [0]*len( prms_file )
dphi_sig = [0]*len( prms_file )

for i in range( len( prms_file ) ) :

	date_str[i] = prms_file[i][0:-1].split(' ')[0]
	date[i]     = datetime( year=int(date_str[i][0:4]), month=int(date_str[i][5:7]),
	                        day=int(date_str[i][8:10]), hour=int(date_str[i][11:13]),
	                        minute=int(date_str[i][14:16]), second=int(date_str[i][17:19]) )
	s[i]        = float(prms_file[i][0:-1].split(' ')[1])
	gA[i]       = float(prms_file[i][0:-1].split(' ')[2])
	gA_sig[i]   = float(prms_file[i][0:-1].split(' ')[3])
	gV[i]       = float(prms_file[i][0:-1].split(' ')[4])
	gV_sig[i]   = float(prms_file[i][0:-1].split(' ')[5])
	dthe[i]     = float(prms_file[i][0:-1].split(' ')[6])
	dthe_sig[i] = float(prms_file[i][0:-1].split(' ')[7])
	dphi[i]     = float(prms_file[i][0:-1].split(' ')[8])
	dphi_sig[i] = float(prms_file[i][0:-1].split(' ')[9])

	#plot_prms()

def plot_prms() :

	f1, axs1 = plt.subplots( 4, 1, sharex=True, squeeze=True )

	axs1[0].set_title( date_str[0].replace(';','/') + ' - ' + date_str[-1].replace(';','/'), fontsize=22 )

	# gA

	#axs1[0].axhline( gA[best_i], c='b', lw='1', label=r'$g_A(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gA[best_i], 4) ) )
	axs1[0].scatter( date, gA, color='k' )
	axs1[0].errorbar( date, gA, yerr=gA_sig, color='k', ls='', elinewidth=2 )
	axs1[0].set_xlim( date[0]-timedelta(minutes=5), date[-1]+timedelta(minutes=5) )
	axs1[0].set_ylabel( r'$g_A$', fontsize=18 )

	# gV

	#axs1[1].axhline( gV[best_i], c='b', lw='1', label=r'$g_V(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(gV[best_i], 4) ) )
	axs1[1].scatter( date, gV, color='k' )
	axs1[1].errorbar( date, gV, yerr=gV_sig, color='k', ls='', elinewidth=2 )
	axs1[1].set_xlim( date[0]-timedelta(minutes=5), date[-1]+timedelta(minutes=5) )
	axs1[1].set_ylabel( r'$g_V$', fontsize=18 )

	# dthe

	#axs1[2].axhline( dthe[best_i], c='b', lw='1', label=r'$\Delta\theta(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dthe[best_i], 4) ) )
	axs1[2].scatter( date, dthe, color='k' )
	axs1[2].errorbar( date, dthe, yerr=dthe_sig, color='k', ls='', elinewidth=2 )
	axs1[2].set_xlim( date[0]-timedelta(minutes=5), date[-1]+timedelta(minutes=5) )
	axs1[2].set_ylabel( r'$\Delta\theta$ $(^\circ)$', fontsize=18 )

	# dphi

	#axs1[3].axhline( dphi[best_i], c='b', lw='1', label=r'$\Delta\phi(s|_{\sigma_{{g_A}_{min}}})=$'+str(round(dphi[best_i], 4) ) )
	axs1[3].scatter( date, dphi, color='k' )
	axs1[3].errorbar( date, dphi, yerr=dphi_sig, color='k', ls='', elinewidth=2 )
	axs1[3].set_xlim( date[0]-timedelta(minutes=5), date[-1]+timedelta(minutes=5) )
	axs1[3].set_ylabel( r'$\Delta\phi$ $(^\circ)$', fontsize=18 )

	axs1[3].set_xlabel( 'Time', fontsize=18 )

	#for tick in axs1[1].yaxis.get_major_ticks():
	#	tick.label.set_fontsize(18)

	#plt.xticks( rotation=45, fontsize=16 )

	#plt.yticks( fontsize=16 )

	plt.subplots_adjust(wspace=0, hspace=0.2)
	plt.tight_layout()

	plot_vars = 'intercal_params_plots_' + date_str[0] + '_' + date_str[-1]

	plt.savefig( os.path.join( 'chi-squared', plot_vars + '.pdf' ), bbox_inches='tight', dpi=40 )
	plt.savefig( os.path.join( 'chi-squared', plot_vars + '.eps' ), bbox_inches='tight', dpi=40 )

	plt.show()
