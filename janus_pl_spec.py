################################################################################
##
## Janus -- GUI Software for Processing Thermal-Ion Measurements from the
##          Wind Spacecraft's Faraday Cups
##
## Copyright (C) 2016 Bennett A. Maruca (bmaruca@udel.edu)
##
## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see http://www.gnu.org/licenses/.
##
################################################################################


################################################################################
## LOAD THE NECESSARY MODULES.
################################################################################

from janus_pl_dat import pl_dat
from datetime import timedelta
from scipy.interpolate import interp1d
from numpy import amax, amin, append, arccos, arctan2, arange, argsort, array, \
                    average, cos, deg2rad, diag, dot, exp, indices, interp, \
                    mean, pi, polyfit, rad2deg, reshape, shape, sign, sin, sum,\
                    sqrt, std, tile, transpose, where, zeros

# Load the modules necessary for signaling the graphical interface.

from PyQt4.QtCore import QObject, SIGNAL, QThread

# Load the dictionary of physical constants.

from janus_const import const

# Load the necessary array modules and mathematical functions.

from numpy import amax, amin, append, arccos, arctan2, arange, argsort, array, \
                    average, cos, deg2rad, diag, dot, exp, indices, interp, \
                    mean, pi, polyfit, rad2deg, reshape, sign, sin, sum, sqrt, \
                    std, tile, transpose, where, zeros

# Load the "pyon" module.

from janus_pyon import plas, series

# Load the modules necessary for copying.

from copy import deepcopy

################################################################################
## DEFINE THE "pl_spec" CLASS.
################################################################################

class pl_spec( ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self,
	              t_strt=None, t_stop=None, elev_cen=None, the_del=None,
	              azim_cen=None, phi_del=None, volt_cen=None, volt_del=None,
	              psd=None                                                   ) :

		self._n_bin  = 14 #TODO Confirm
		self._n_the  = 5
		self._n_phi  = 5  #TODO Confirm this
		self._t_strt = t_strt
		self._t_stop = t_stop
		self._rot    = t_stop - t_strt

		if ( azim_cen is None ) :
			azim_cen = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			azim_cen = self.adjust(azim_cen)

		if ( phi_del is None ) :
			phi_del  = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			phi_del = self.adjust(phi_del)

		if ( elev_cen is None ) :
			elev_cen = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			elev_cen = self.adjust(elev_cen)


		if ( the_del is None ) :
			the_del  = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			the_del = self.adjust(the_del)

		if ( volt_cen is None ) :
			volt_cen = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			volt_cen = self.adjust(volt_cen)	

		if ( volt_del is None ) :
			volt_del = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			volt_del = self.adjust(volt_del)	

		if ( psd is None ) :
			psd      = [ [ [ None for b in range( self._n_bin ) ]
			                      for p in range( self._n_phi ) ]
			                      for t in range( self._n_the ) ]
		else:
			psd = self.adjust(psd)

		self.arr = [[[ pl_dat( spec=self,
		                       azim_cen = azim_cen[t][p][b],
		                       phi_del  = phi_del[t][p][b],
		                       elev_cen = elev_cen[t][p][b],
		                       the_del  = the_del[t][p][b],
		                       volt_cen = volt_cen[t][p][b], 
		                       volt_del = volt_del[t][p][b], 
		                       t_strt   = self._t_strt,
		                       t_stop   = self._t_stop,
		                       psd      = psd[t][p][b]       ) 
		                              for b in range( self._n_bin ) ]
		                              for p in range( self._n_phi ) ]
		                              for t in range( self._n_the ) ]

		# Validate the data in the spectrum.

#		self.validate( )

		# List of shape [[[n_bin] n_the] n_phi] where
		# 'bin' is the voltage sweep number,
		# 'the' specifies theta of the look direction, and
		# 'phi' specifies phi of the look direction

		# Initialize the variables that will contain the settings,
		# data selections, and results from all analyses.

		self.rset_spec( var_spec    = False,
		                var_mom_win = True,
		                var_mom_sel = True,
		                var_mom_res = True,
		                var_nln_ion = True,
		                var_nln_set = True,
		                var_nln_gss = True, 
		                var_nln_sel = True,
		                var_nln_res = True,
		                var_dsp     = True,
		                var_dyn     = True  )
		
		# Initialize the value of the indicator variable of whether the
		# automatic analysis should be aborted.

		# Note.  This variable provides the user with a means of
		#        prematurely stopping the automatic analysis of a series
		#        of spectra by the function "self.auto_run".  When that
		#        function is called, this indicator's value is set to
		#        "False".  After the function processes a given
		#        spectrum, it procedes on to the next one only if this
		#        indicator is still "False".

		self.stop_auto_run = False

	def adjust(self, matrix) :

	# This function 'adjust':
	# 1. takes data whose phi values and voltages are reversed and returns
	# them to their proper order (it does not affect theta values) and
	# 2. breaks 25 directions into 5 sets of 5 pairs of (theta, phi)

		new_matrix = [ [ [ None for b in range( self._n_bin ) ]
		                        for p in range( self._n_phi ) ]
		                        for t in range( self._n_the ) ]
		for t in range( self._n_the ):
			for p in range( self._n_phi ):
				for b in range( self._n_bin ):
					new_matrix[-t-1][p][b] = \
					    matrix[-(5*t+p)-1][-b-1]
		return new_matrix

	def __getitem__(self, key ) : #TODO not yet finished

		if ( key == 'n_the' ) :
			return self._n_the
		elif ( key == 'n_phi' ) :
			return self._n_phi
		elif ( key == 'n_bin' ) :
			return self._n_bin
		elif ( key == 'n_dir' ) :
			return ( self._n_the * self._n_phi )
		elif ( key == 'time' ) :
			return [ self.arr[0][p][0]['time']
					for p in range(self._n_phi)]
		elif ( key == 'elev_cen' ) :
			return [ self.arr[t][0][0]['elev_cen']
					for t in range(self._n_the)]
		elif ( key == 'the_del' ) :
			return [ self.arr[t][0][0]['the_del']
					for t in range(self._n_the)]
		elif ( key == 'azim_cen' ) :
			return [ self.arr[0][p][0]['azim_cen']
					for p in range(self._n_phi)]
		elif ( key == 'phi_del' ) :
			return [ self.arr[0][p][0]['phi_del']
					for p in range(self._n_phi)]
		elif ( key == 'volt_strt' ) :
			return  [ self.arr[0][0][b]['volt_strt'] 
					for b in range(self._n_bin)]
		elif ( key == 'volt_stop' ) :
			return  [ self.arr[0][0][b]['volt_stop'] 
					for b in range(self._n_bin)]
		elif ( key == 'volt_cen' ) :
			return  [ self.arr[0][0][b]['volt_cen'] 
					for b in range(self._n_bin)]
		elif ( key == 'volt_del' ) :
			return  [ self.arr[0][0][b]['volt_del'] 
					for b in range(self._n_bin)]
		elif ( key == 'vel_strt' ) :
			return  [ self.arr[0][0][b]['vel_strt'] 
					for b in range(self._n_bin)]
		elif ( key == 'vel_stop' ) :
			return  [ self.arr[0][0][b]['vel_stop'] 
					for b in range(self._n_bin)]
		elif ( key == 'vel_cen' ) :
			return  [ self.arr[0][0][b]['vel_cen'] 
					for b in range(self._n_bin)]
		elif ( key == 'vel_del' ) :
			return  [ self.arr[0][0][b]['vel_del'] 
					for b in range(self._n_bin)]
		elif ( key == 'psd' ) :
			return [ [ [ self.arr[t][p][b]['psd']
			             for b in range( self._n_bin ) ]
			             for t in range( self._n_phi ) ]
			             for p in range( self._n_the ) ]
		elif ( key == 'psd_valid' ) :
			return [ [ [ self.arr[t][p][b]['psd']
			             for b in range( self._n_bin ) ]
			             for t in range( self._n_phi ) ]
			             for p in range( self._n_the ) ]
		elif ( key == 'psd_flat' ) :
			return [ self.arr[t][p][b]['psd']
			         for b in range( self._n_bin )
			         for t in range( self._n_phi )
			         for p in range( self._n_the ) ]
		elif ( key == 'rot' ) :
			return self._rot
		else :
			raise KeyError( 'Invalid key for "pl_spec".' )

	def __setitem__( self, key, val ) :		

		raise KeyError( 'Reassignment not permitted except through'
		                                    + ' "set_*" functions.' )

        #-----------------------------------------------------------------------
	# DEFINE THE FUNCTION TO ASSIGN THE MAGNETIC FIELD TO EACH DATUM. 
	#-----------------------------------------------------------------------

	def set_mag( self, mfi_t, mfi_b_x, mfi_b_y, mfi_b_z ) :

		mfi_s = [ ( t - mfi_t[0] ).total_seconds( ) for t in mfi_t ]

		fnc_b_x = interp1d( mfi_s, mfi_b_x )
		fnc_b_y = interp1d( mfi_s, mfi_b_y )
		fnc_b_z = interp1d( mfi_s, mfi_b_z )

		try :

			for t in range( self['n_the'] ) :

				for p in range( self['n_phi'] ) :

					for b in range( self['n_bin'] ) :

						s = ( self.arr[t][p][b]['time']
                                                  - mfi_t[0] ).total_seconds( )

						b_x = fnc_b_x( s )
						b_y = fnc_b_y( s )
						b_z = fnc_b_z( s )

						self.arr[t][p][b].set_mag( (
						            b_x, b_y, b_z ) )
		except :

			avg_b_x = sum( mfi_b_x ) / float( len( mfi_b_x ) )
			avg_b_y = sum( mfi_b_y ) / float( len( mfi_b_y ) )
			avg_b_z = sum( mfi_b_z ) / float( len( mfi_b_z ) )

			for t in range( self['n_the'] ) :

                                for p in range( self['n_phi'] ) :

                                        for b in range( self['n_bin'] ) :

                                                self.arr[t][p][b].set_mag( (
						avg_b_x, avg_b_y, avg_b_z ) )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RUNNING THE MOMENTS ANALYSIS ON THIS SPECTRUM.
	#-----------------------------------------------------------------------

	def anls_mom( self ) :

		# Re-initialize and the output of the moments analysis.

		self.rset_spec( var_mom_res=True )

		# If the point-selection arrays have not been populated, run
		# the automatic point selection.

		if ( ( self.mom_sel_dir is None ) or
		     ( self.mom_sel_bin is None )    ) :

			self.auto_mom_sel( )

		# If any of the following conditions are met, emit a
		# signal that indicates that the results of the moments
		#analysis have changed, and then abort.
		#   -- No (valid) ion spectrum has been requested.
		#   -- Insufficient data have been selected.

		if ( ( self.arr is None  ) or
		     ( self.mom_n_sel_dir < self.mom_min_sel_dir ) or
		     ( self.mom_n_sel_dir > self['n_dir']        )    ) :
			# Note:  n_dir is set to the total number of
			#        (theta, phi) directions
			self.emit( SIGNAL('janus_mesg'),
			                  'core', 'norun', 'mom' )

			self.emit( SIGNAL('janus_chng_mom_res') )

			return

		# Extract the "t"- and "p"-indices of each selected
		# pointing direction.

		( tk_t, tk_p ) = where( self.mom_sel_dir )

		# Initialize the "eta_*" arrays.

		# Note.  Only some of these arrays will ultimately be
		# saved.

		# Note.  The arrays "eta_?" define the density,
		#        inflow speed, inflow velocity, thermal speed,
		#        and temperature derived for each of the
		#        analyzed look directions.

		n_eta = self.mom_n_sel_dir

		eta_dlk = tile( 0., [ n_eta, 3 ] )     # Cartesian look
	                                               # direction
		eta_n     = tile( 0., n_eta )          # number density
		eta_v     = tile( 0., n_eta )          # inflow speed
		eta_v_vec = tile( 0., [ n_eta, 3 ] )   # inflow velocity
		eta_w     = tile( 0., n_eta )          # thermal speed
		eta_t     = tile( 0., n_eta )          # temperature

		# For each of the selected look directions, identify the
		# selected data therefrom and calculate the estimator of
		# the projected inflow speed along that direction.

		for k in range( n_eta ) :

			# Extract the "t"- and "p"-values for this
			# direction.

			t = tk_t[k]
			p = tk_p[k]

			# Calculate the look direction using "t"- and
			# "p"-values

			eta_dlk[k] = self.arr[t][p][0]['dir']

			# Extract the "b" values of the selected data
			# from this look direction.

                        b = [i for i, x in enumerate(
			                    self.mom_sel_bin[t][p])
                                                          if x==True   ]

			# Define the variables for this calculation

			this_win = self.arr[t][p]

			f_u     = [ this_win[i]['psd']     for i in b ]

			u       = [ this_win[i]['vel_cen'] for i in b ]

			u_vec   = [ this_win[i]['dir']     for i in b ]

			d_u     = [ this_win[i]['vel_del'] for i in b ]

			theta   = [ deg2rad( this_win[i]['the_cen'] )
			                                    for i in b ]

			d_theta = [ deg2rad( this_win[i]['the_del'] )
			                                    for i in b ]

			d_phi   = [ deg2rad( this_win[i]['phi_del'] )
			                                    for i in b ]

			d_omega = [ sin( theta[i] ) * d_theta[i] *
			            d_phi[i] for i in range(len(b))         ]

			# Compute the number density for this spectrum.

			eta_n[k] = sum( [ f_u[i] * u[i]**2 * d_u[i] *
			                  d_omega[i] for i in range(len(b)) ]     )

			# Compute the bulk velocity.

			eta_v_vec[k] = [ sum([f_u[i] * u_vec[i][j] *
			                      u[i]**2 * d_u[i] *
			                      d_omega[i] for i in range(len(b)) ])/
			                 eta_n[k] for j in range(3)    ]

			# Compute the bulk speed.

			eta_v[k] = sqrt( sum( [ eta_v_vec[k][j]**2
			                         for j in range(3) ] ) )

			# Compute the thermal speed.

			eta_w[k] = sqrt( ( sum( [ f_u[i] * u[i]**4 *
			                          d_u[i] * d_omega[i]
			                            for i in range(len(b)) ]     )/
			                 eta_n[k] - eta_v[k]**2 ) / 3. )

			# Compute the effective temperature.

			eta_t = ( 1.E-3 / const['k_b'] ) * \
		        	const['m_p'] * ( ( 1.E3  * eta_w )**2 )

		# Calculate a net estimators of the number density and
		# thermal speed.

		mom_n = mean( eta_n )

		mom_w = mean( eta_w )

		mom_v_vec = [ mean( eta_v_vec[j] ) for j in range(3) ]

		# Save the results of the moments analysis in a plas object.

		self.mom_res = plas( )

		self.mom_res['v0_vec'] = mom_v_vec

		self.mom_res.add_spec( name='Proton', sym='p', m=1., q=1. )

		self.mom_res.add_pop( 'p',
		                      drift=False, aniso=False,
		                      name='Core', sym='c',
		                      n=mom_n,     w=mom_w          )

		# Calculate the expected PSDs based on the results of the
		# (linear) moments analysis.

		# TODO

		return self.mom_res

###########################################################################################################

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR AUTOMATIC DATA SELECTION FOR THE MOMENTS ANLS.
	#-----------------------------------------------------------------------

	def auto_mom_sel( self ) :

		# Re-initialize the data-selection variables for the moments
		# analysis.

		self.rset_spec( var_mom_sel=True )

		# If no spectrum has been loaded, abort.

		if ( self.arr is None ) :

			return

		# Initially, deselect all look directions and bins.

		self.mom_sel_dir = [ [ False 
		                       for p in range( self['n_phi'] ) ]
		                       for t in range( self['n_the'] ) ]

		self.mom_sel_bin = [ [ [ False 
		                         for b in range( self['n_bin'] ) ]
                                         for p in range( self['n_phi'] ) ]
		                         for t in range( self['n_the'] ) ]

		# If the "mom_win_???" variables are invalid, abort.

		if ( ( self.mom_win_dir is None ) or
		     ( self.mom_win_bin is None )    ) :

			self.vldt_mom_sel( emit_all=True )

			self.anls_mom( )

			return

		# Find the maximum psd window (of "self.mom_win_bin" bins)
		# for each direction
		dir_max_ind  = [ [ self.find_max_psd( t, p,
		                             win=self.mom_win_bin              )
		                        for p in range(self['n_phi'] ) ]
		                        for t in range(self['n_the'] ) ]

		dir_max_psd = [ [ self.calc_tot_psd( t, p,
		                             dir_max_ind[t][p],
		                             win=self.mom_win_bin              )
		                        for p in range(self['n_phi'] ) ]
		                        for t in range(self['n_the'] ) ]

		# Compute the "self.mom_win_dir" window indices with the highest
		# total PSD

		win_max_ind = []

		for t in range( self['n_the'] ) :

			for p in range( self['n_phi'] ) :

				n_big = 0

				for tp in range( self['n_dir'] ) :

					if ( dir_max_psd[t][p] < dir_max_psd[ tp//self['n_the'] ][ tp % self['n_the'] ] ) :
						n_big += 1

				if ( n_big < self.mom_win_dir ) :
					win_max_ind += [[t, p]]

				if len( win_max_ind ) == self.mom_win_dir :

					break

			if len( win_max_ind ) == self.mom_win_dir :

				break

		# Populate "self.mom_sel_bin" and "self.mom_sel_dir"
		# appropriately.

		for win in win_max_ind :

			t = win[0]
			p = win[1]

			self.mom_sel_dir[t][p] = True

			# Select the bins in this look direction's
			# maximal window

			for b in range( dir_max_ind[t][p],
			                dir_max_ind[t][p] + self.mom_win_bin ) :
				self.mom_sel_bin[t][p][b] = True

                # Validate the new data selection (which includes populating
		# the "self.mom_sel_dir" array).

		self.vldt_mom_sel( emit_all=True )

		# Run the moments analysis (and then, if the non-linear analysis
		# is set to be dynamically updated, run that analysis as well).

		self.anls_mom( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION TO FIND THE INDEX OF WINDOW WITH MAXIMUM PSD
	#-----------------------------------------------------------------------

	def find_max_psd( self, t, p, win=1 ) :

		# Validate the theta and phi indices.

		if ( ( t < 0 ) or ( t >= self['n_the'] ) ) :
			raise ValueError( 'Theta index out of range.' )

		if ( ( p < 0 ) or ( p >= self['n_phi'] ) ) :
			raise ValueError( 'Phi index out of range.' )

		# Validate the window size.

		if ( ( win < 1 ) or ( win > self['n_bin'] ) ) :
			raise ValueError( 'Window out of range.' ) 

		# Search the specified direction for the "win"-bin range that
		# contains the maximum total psd.

		b_max   = 0
		psd_max = 0.

		for b in range( 0, self['n_bin'] - win + 1 ) :

			psd = sum( [ self.arr[t][p][b+w]['psd']
			                               for w in range( win ) ] )

			if ( psd > psd_max ) :
				b_max    = b
				psd_max = psd

		# Return the location of the window with the maximum psd.

		return b_max

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR CALCULATING TOTAL CURRENT IN A GIVEN WINDOW.
	#-----------------------------------------------------------------------

	def calc_tot_psd( self, t, p, b, win=1 ) :

		# Validate the theta, phi, and bin indices.

		if ( ( t < 0 ) or ( t >= self['n_the'] ) ) :
			raise ValueError( 'Theta index out of range.' )

		if ( ( p < 0 ) or ( p >= self['n_phi'] ) ) :
			raise ValueError( 'Phi index out of range.' )

		if ( ( b < 0 ) or ( b >= self['n_bin'] ) ) :
			raise ValueError( 'Velocity index out of range.' )

		# Validate the window size.

		if ( ( win < 1 ) or ( b + win > self['n_bin'] ) ) :
			raise ValueError( 'Window out of range.' ) 

		# Return the total (valid) current in the specified window.

		return sum( [ self.arr[t][p][b+w]['psd']
		                                       for w in range( win ) ] )

	def vldt_mom_sel( self, emit_all=False ) :

		# Note.  This function ensures that the two "self.mom_sel_???"
		#        arrays are mutually consistent.  For each set of "c"-
		#        and "d"-values, "self.mom_sel_dir[c,d]" can only be
		#        "True" if at least "self.min_sel_bin" of the elements
		#        in "self.mom_sel_bin[c,d,:]" are "True".  However, if
		#        fewer than "self.mom_min_sel_dir" sets of "c"- and
		#        "d"-values satisfy this criterion, all elements of
		#        "self.mom_sel_dir" are given the value "False".		
		#
		#        Additionally, this functions serves to update the
		#        "self.mom_n_sel_???" counters.


		# Save the initial selection of pointing directions.

		old_mom_sel_dir = deepcopy( self.mom_sel_dir )

		# Update the counter "self.mom_n_sel_bin" (i.e., the number of
		# selected data in each pointing direction).

		self.mom_n_sel_bin = [ [ sum( self.mom_sel_bin[t][p] )
		                       for p in range( self['n_phi'] ) ]
		                       for t in range( self['n_the'] ) ]

		# Create a new selection of pointing directions based on the
		# data selection, and then update the counter
		# "self.mom_n_sel_dir".

		self.mom_sel_dir = [ [
		              self.mom_n_sel_bin[t][p] >= self.mom_min_sel_bin
		                       for p in range( self['n_phi'] ) ]
		                       for t in range( self['n_the'] ) ]

		# Determine the total number of selected pointing directions; if
		# this number is less than the minimum "self.mom_min_sel_dir",
		# deselect all pointing directions.

		self.mom_n_sel_dir = \
		               sum( [ sum( sub ) for sub in self.mom_sel_dir ] )

		if ( self.mom_n_sel_dir < self.mom_min_sel_dir ) :

			self.mom_sel_dir = [ [ False
			               for p in range( self['n_dir'] ) ]
			               for t in range( self['n_cup'] ) ]

			self.mom_n_sel_dir = 0

		"""
		# Emit (if necessary) the appropriate update signal(s).

		if ( emit_all ) :

			self.emit( SIGNAL('janus_chng_mom_sel_all') )

		else :

			# Identify differences between the new and old versions
			# of "self.mom_sel_dir".  For each pointing direction
			# whose selection status for the moments analysis has
			# changed, emit a signal indicating this.

			for t in range( self['n_cup'] ) :

				for p in range( self['n_dir'] ) :

					if ( self.mom_sel_dir[t][p]
					            != old_mom_sel_dir[t][p] ) :

						self.emit( SIGNAL(
						      'janus_chng_mom_sel_dir'),
						      t, p )
		"""

	#-----------------------------------------------------------------------
	# RESET THE DATA AND ANALYSIS VARIABLES.
	#-----------------------------------------------------------------------

	def rset_spec( self,
	              var_spec    = False,
	              var_mom_win = False,
		      var_mom_sel = False,
	              var_mom_res = False,
	              var_nln_ion = False,
	              var_nln_set = False,
		      var_nln_gss = False,
	              var_nln_sel = False,
		      var_nln_res = False,
	              var_dsp     = False,
	              var_dyn     = False  ) :

		# If requested, (re-)initialize the varaibles for the Wind/PESA
		# data associated with this spectrum.

		if ( var_spec ) :

			self._n_bin  = 14 #TODO Confirm
			self._n_the  = 5
			self._n_phi  = 5  #TODO Confirm this
			self._t_strt = None
			self._t_stop = None
			self._rot    = None
			self.arr     = None

		# If requested, (re-)initialize the varaibles for the windows
		# associated with automatic data selection for the PL moments
		# analysis.

		if ( var_mom_win ) :

			self.mom_win_dir = 7
			self.mom_win_bin = 7

		# If requested, (re-)initialize the variables associated with
		# the data seleciton for the PL moments analysis.

		if ( var_mom_sel ) :

			self.mom_min_sel_dir =  5
			self.mom_min_sel_bin =  3

			self.mom_max_sel_dir = 25

			self.mom_sel_dir     = None
			self.mom_sel_bin     = None

		# If requested, (re-)initialize and store the variables
		# associated with the results of the PL moments analysis.

		if ( var_mom_res ) :

			self.mom_res  = None

		# If requested, (re-)initialize the variables associated with
		# the ion species and populations for the non-linear analysis.

		# Note.  This includes both the "self.nln_spc_?" and
		#        "self.nln_pop_?" arrays.  These are done together since
		#        they are so interconnected by the "self.nln_plas"
		#        object, which is also handled here.

		"""
		if ( var_nln_ion ) :

			self.nln_n_spc = 4
			self.nln_n_pop = 5

			self.nln_plas = plas( enforce=True )

			self.nln_pop_use = tile( False, self.nln_n_pop )
			self.nln_pop_vld = tile( False, self.nln_n_pop )

			for s in range ( self.nln_n_spc ) :

				if ( s == 0 ) :
					self.nln_plas.add_spec( name='Proton',
					                   sym='p', m=1., q=1. )
				elif ( s == 1 ) :
					self.nln_plas.add_spec( name='Alpha' ,
					                   sym='a', m=4., q=2. )
				else :
					self.nln_plas.add_spec( )

			for p in range ( self.nln_n_pop ) :

				if ( p == 0 ) :
					self.nln_pop_use[p] = True
					self.nln_pop_vld[p] = True
					self.nln_plas.add_pop(
					        'p', name='Core', sym='c',
					        drift=False, aniso=True    )
				elif ( p == 1 ) :
					self.nln_pop_use[p] = False
					self.nln_pop_vld[p] = True
					self.nln_plas.add_pop(
					        'p', name='Beam', sym='b',
					        drift=True , aniso=False   )
				elif ( p == 2 ) :
					self.nln_pop_use[p] = True
					self.nln_pop_vld[p] = True
					self.nln_plas.add_pop(
					        'a', name='Core', sym='c',
					        drift=True , aniso=True    )
				elif ( p == 3 ) :
					self.nln_pop_use[p] = False
					self.nln_pop_vld[p] = True
					self.nln_plas.add_pop(
					        'a', name='Beam', sym='b',
					        drift=True , aniso=False   )
				else :
					self.nln_pop_use[p] = False
					self.nln_pop_vld[p] = False
					self.nln_plas.add_pop( None )

		# If requested, (re-)initialize the variables associated with
		# the settings for the automatic initial guess generation and
		# the automatic point selection.

		if ( var_nln_set ) :

			self.nln_set_gss_n   = tile( None , self.nln_n_pop )
			self.nln_set_gss_d   = tile( None , self.nln_n_pop )
			self.nln_set_gss_w   = tile( None , self.nln_n_pop )
			self.nln_set_gss_vld = tile( False, self.nln_n_pop )

			self.nln_set_sel_a   = tile( None , self.nln_n_pop )
			self.nln_set_sel_b   = tile( None , self.nln_n_pop )
			self.nln_set_sel_vld = tile( False, self.nln_n_pop )

			self.nln_set_gss_n[0]   =  1.00
			self.nln_set_gss_n[1]   =  0.20
			self.nln_set_gss_n[2]   =  0.02
			self.nln_set_gss_n[3]   =  0.01

			self.nln_set_gss_d[1]   =  0.03
			self.nln_set_gss_d[2]   =  0.01
			self.nln_set_gss_d[3]   =  0.05

			self.nln_set_gss_w[0]   =  1.00
			self.nln_set_gss_w[1]   =  1.25
			self.nln_set_gss_w[2]   =  1.00
			self.nln_set_gss_w[3]   =  1.25

			self.nln_set_gss_vld[0] = True
			self.nln_set_gss_vld[1] = True
			self.nln_set_gss_vld[2] = True
			self.nln_set_gss_vld[3] = True

			self.nln_set_sel_a[0]   = -3.00
			self.nln_set_sel_a[1]   = -3.00
			self.nln_set_sel_a[2]   = -3.00
			self.nln_set_sel_a[3]   = -3.00
			self.nln_set_sel_a[4]   = -3.00

			self.nln_set_sel_b[0]   =  3.00
			self.nln_set_sel_b[1]   =  3.00
			self.nln_set_sel_b[2]   =  3.00
			self.nln_set_sel_b[3]   =  3.00
			self.nln_set_sel_b[4]   =  3.00

			self.nln_set_sel_vld[0] = True
			self.nln_set_sel_vld[1] = True
			self.nln_set_sel_vld[2] = True
			self.nln_set_sel_vld[3] = True
			self.nln_set_sel_vld[4] = True

		# If requested, (re-)initialize the variables associated with
		# the initial guesses for the non-linear analysis.

		if ( var_nln_gss ) :

			for p in range( self.nln_n_pop ) :
				self.nln_plas.arr_pop[p]['n']     = None
				self.nln_plas.arr_pop[p]['dv']    = None
				self.nln_plas.arr_pop[p]['w']     = None
				self.nln_plas.arr_pop[p]['w_per'] = None
				self.nln_plas.arr_pop[p]['w_par'] = None

			self.nln_gss_vld = tile( False, self.nln_n_pop )

			self.nln_gss_pop      = [ ]
			self.nln_gss_prm      = [ ]

			self.nln_gss_curr_tot = None
			self.nln_gss_curr_ion = None

		# If requested, (re-)initialize the variables associated with
		# the data selection for the non-linear analysis.

		if ( var_nln_sel ) :

			self.nln_sel     = None

			self.nln_n_sel   = 0
			self.nln_min_sel = 30

		# If requested, (re-)initialize the variables associated with
		# the results of the non-linear analysis.

		if ( var_nln_res ) :

			self.nln_res_runtime  = 0.

			self.nln_res_plas     = None

			self.nln_res_sel      = None

			self.nln_res_curr_tot = None
			self.nln_res_curr_ion = None

		"""
		# If requested, (re-)initialize the variables which indicate of
		# the analyses have their results displayed in widgets which
		# support output from multiple analyses.


		if ( var_dsp ) :

			self.dsp = 'mom'

		# If requested, (re-)initialize the variables which indicate
		# which analyses are updated automatically when a change is
		# made to their settings.

		if ( var_dyn ) :

			self.dyn_mom = True
			self.dyn_gss = True
			self.dyn_sel = True
			self.dyn_nln = False
