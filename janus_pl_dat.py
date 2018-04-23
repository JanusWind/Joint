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

from math import sqrt, acos, pi
from numpy import interp, sin, cos, deg2rad, exp, array
from scipy.special import erf
from datetime import datetime, timedelta

from janus_const import const
from janus_helper import calc_arr_norm, calc_arr_dot


################################################################################
## DEFINE THE CLASS FOR DATUM
################################################################################

class pl_dat( ) :

	def __init__( self, spec=None,
	              t_strt=None, t_stop=None, azim_cen=None, phi_del=None,
	              elev_cen=None, the_del=None, volt_cen=None, volt_del=None,
	              psd=None, valid=False, mom_sel = False ) :

		self._spec      = spec
		self._azim_cen  = azim_cen
		self._phi_del   = phi_del
		self._elev_cen  = elev_cen
		self._the_del   = the_del
		self._volt_cen  = volt_cen
		self._volt_del  = volt_del
		self._psd       = psd
		self._valid     = valid
		self.mom_sel    = mom_sel

		self._time = ( datetime( 1970, 1, 1 ) + 
		               timedelta( seconds = (t_strt +
		               ( t_stop - t_strt ) / 360. * azim_cen ) ) )

		#Note: The voltage is recorded from high to low voltage
		self._volt_strt = ( self._volt_cen + ( self._volt_del / 2. ) )
		self._volt_stop = ( self._volt_cen - ( self._volt_del / 2. ) )

		self._vel_strt  = 1E-3*( sqrt(2.0*const['q_p']*
		                         self['volt_strt']/const['m_p'])    )
		self._vel_stop  = 1E-3*( sqrt((2.0*const['q_p']*
		                         self['volt_stop']/const['m_p']))   )

		self._vel_cen   = 1E-3*( sqrt(2.0*const['q_p']*
		                         self._volt_cen/const['m_p'])       )

		self._vel_del   = (  self['vel_strt']-self['vel_stop']      )

		# TODO It is currently assumed that the given values of theta
		#      and phi are the proper look directions. Need to confirm.

		# TODO: Confirm these formulae

		self._the       =( 90 + self._elev_cen ) * pi/180 # ( 90 -
		self._phi       =( 180 + self._azim_cen) * pi/180

		self._dir_x     = - sin( self._the ) * cos( self._phi )
		self._dir_y     = - sin( self._the ) * sin( self._phi )
		self._dir_z     = - cos( self._the )
		#/TODO

		self._norm_b_x  = None
		self._norm_b_y  = None
		self._norm_b_z  = None

		if ( ( self._time     is None ) or
		     ( self._azim_cen is None ) or ( self._phi_del  is None ) or
                     ( self._elev_cen is None ) or ( self._the_del  is None ) or
		     ( self._volt_cen is None ) or ( self._volt_del is None ) or
		     ( self._psd      is None )                              ) :
			self._valid = False
		else :
			self._valid = True

		# Define the variables for use in the moments analysis

		self._mom0  = self['psd']*self['vel_cen']**2*\
		              sin( self['the_cen'] )*self['vel_del']*\
		              deg2rad( self['the_del'] )*\
		              deg2rad( self['phi_del'] )
		self._mom1x = self._mom0*self['vel_cen']*self._dir_x
		self._mom1y = self._mom0*self['vel_cen']*self._dir_y
		self._mom1z = self._mom0*self['vel_cen']*self._dir_z
		self._mom2  = self._mom0*self['vel_cen']**2

	# ----------------------------------------------------------------------
	# DEFINE THE KEYS FOR THIS CLASS.
	# ----------------------------------------------------------------------

	def __getitem__( self, key ) :
#
#               return self.__dict__['_'+key]
#
		if ( key == 'spec' ) :
			return self._spec
		elif ( key == 'valid' ) :
			return self._valid
		elif ( key == 't_strt' ) :
			return self._t_strt
		elif ( key == 't_stop' ) :
			return self._t_stop
                elif ( key == 'time' ) :
                        return self._time
		elif ( key == 'volt_cen' ) :
			return self._volt_cen
		elif ( key == 'volt_del' ) :
			return self._volt_del
		elif ( key == 'volt_strt' ) :
			return self._volt_strt
		elif ( key == 'volt_stop' ) :
			return self._volt_stop
		elif ( key == 'vel_strt' ) :
			return self._vel_strt
		elif ( key == 'vel_stop' ) :
			return self._vel_stop
		elif ( key == 'vel_cen' ) :
			return self._vel_cen
		elif ( key == 'vel_del' ) :
			return self._vel_del
		elif ( key == 'psd' ) :
			return self._psd
		elif ( key == 'psd_valid' ) :
			if ( self['valid'] ) :
				return self['psd']
			else :
				return 0.
		elif ( key == 'azim_cen' ) :
			return self._azim_cen
		elif ( key == 'elev_cen' ) :
			return self._elev_cen
		elif ( key == 'the_cen' ) :
			return self._the
		elif ( key == 'the_del' ) :
			return self._the_del
		elif ( key == 'phi_cen' ) :
			return self._phi
		elif ( key == 'phi_del' ) :
			return self._phi_del
		elif ( key == 'dir_x' ) :
			return self._dir_x
		elif ( key == 'dir_y' ) :
			return self._dir_y
		elif ( key == 'dir_z' ) :
			return self._dir_z
		elif ( key == 'dir' ) :
			return ( self._dir_x, self._dir_y, self._dir_z )
		elif ( key == 'norm_b_x' ) :
			return self._norm_b_x
		elif ( key == 'norm_b_y' ) :
			return self._norm_b_y
		elif ( key == 'norm_b_z' ) :
			return self._norm_b_z
		elif ( key == 'norm_b' ) :
			return ( self._norm_b_x,self._norm_b_y,self._norm_b_z )
		elif ( key == 'mom_sel' ) :
			return ( self.mom_sel )
		elif ( key == 'mom0' ) :
			return ( self._mom0 )
		elif ( key == 'mom1x' ) :
			return ( self._mom1x )
		elif ( key == 'mom1y' ) :
			return ( self._mom1y )
		elif ( key == 'mom1z' ) :
			return ( self._mom1z )
		elif ( key == 'mom2' ) :
			return ( self._mom2 )
		elif ( key == 'mom0_sel' ) :
			if self.mom_sel :
				return self['mom0']
			else :
				return 0
		elif ( key == 'mom1x_sel' ) :
			if self.mom_sel :
				return self['mom1x']
			else :
				return 0
		elif ( key == 'mom1y_sel' ) :
			if self.mom_sel :
				return self['mom1y']
			else :
				return 0
		elif ( key == 'mom1z_sel' ) :
			if self.mom_sel :
				return self['mom1z']
			else :
				return 0
		elif ( key == 'mom1_sel' ) :
			if self.mom_sel :
				return [ self['mom1x'], self['mom1y'],
				                self['mom1z']          ]
			else :
				return [ 0, 0, 0 ]
		elif ( key == 'mom2_sel' ) :
			if self.mom_sel :
				return self['mom2']
			else :
				return 0
		else :
			raise KeyError( 'Invalid key for "pl_dat".' )

	def __setitem__( self, key, val ) :

		raise KeyError('Reassignment not allowed after initialization.')


	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR SETTING THE MAGNETIC FIELD DIRECTION.
	#-----------------------------------------------------------------------

	def set_mag( self, b_vec ) :

		# Normalize the magnetic-field vector.

		norm_b = calc_arr_norm( b_vec )

		# Store the components of the normalized magnetic-field vector.

		self._norm_b_x = norm_b[0]
		self._norm_b_y = norm_b[1]
		self._norm_b_z = norm_b[2]

		# Compute perpendicular and parallel velocities.

		self._u_par   = self['vel_cen'] * calc_arr_dot( self['norm_b'],
		                                                self['dir']    )
		self._u_par_x = self._u_par * self._norm_b_x
		self._u_par_y = self._u_par * self._norm_b_y
		self._u_par_z = self._u_par * self._norm_b_z

		self._u_per_x = self['vel_cen'] * self['dir_x'] - self._u_par_x
		self._u_per_y = self['vel_cen'] * self['dir_y'] - self._u_par_y
		self._u_per_z = self['vel_cen'] * self['dir_z'] - self._u_par_z
		self._u_per   = sqrt( self._u_per_x**2 + self._u_per_y**2 +
		                                         self._u_per_z**2  )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR SETTING THE MOMENTS SELECTION BOOLEAN.
	#-----------------------------------------------------------------------

	def set_mom_sel( self, sel ) :

		self.mom_sel = sel
