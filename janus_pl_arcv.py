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

# Load the necessary modules for signaling the graphical interface.

from PyQt4.QtCore import SIGNAL

# Load the modules necessary handling dates and times.

from datetime import datetime, timedelta

from janus_time import calc_time_str, calc_time_val, calc_time_epc

from janus_pl_spec import pl_spec

# Load the necessary "numpy" array modules.

from numpy import argsort, array, where, zeros

from operator import attrgetter		  

# Load the modules necessary for file I/O (including FTP).

from scipy.io import readsav

import os.path

from glob import glob

from ftplib import FTP


################################################################################
## DEFINE THE CLASS pl_tag TO HAVE SPECTRA FOR PARTICULAR TIME STAMP
################################################################################

class pl_tag() :
	
	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, c=None, s=None, epoch=None ) :

		self.c     = c
		self.s     = s
		self.epoch = epoch


################################################################################
## DEFINE THE "pl_arcv" CLASS FOR ACCESSING THE ARCHIVE OF Wind/PL SPECTRA.
################################################################################

class pl_arcv( object ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core=None, buf=None, tol=None,
	                    path=None, verbose=None                   ) :

		# Save the arguments for later use.

		# Note.  The "n_file_max" argument is handled at the end of this
		#        function with a call of "chng_n_file_max".

		self.core       = core

		self.buf        = float( buf )      if ( buf is not None )\
		                                    else 3600.

		self.tol        = float( tol )      if ( tol is not None )\
		                                    else 3600.

		self.path       = str( path )       if ( path  is not None )\
		                                    else os.path.join( 
		                                    os.path.dirname( __file__ ), 
		                                    'data', 'pl' )

		self.verbose    = bool( verbose )   if ( verbose is not None )\
		                                    else True

		# Validate the values of parameters.

		if ( self.buf < 0 ) :
			raise ValueError( 'Time buffer cannot be negative.'    )

		# Initialize arrays of date and times loaded.

		self.arr_cdf  = [ ]
		self.arr_date = [ ]

		self.arr_tag  = [ ]

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR LOADING (AND RETURNING) AN ION SPECTRUM.
	#-----------------------------------------------------------------------

	def load_spec( self, time, get_prev=False, get_next=False ) :	

		# If both "get_????" keywords are "True", abort.

		if ( ( get_prev ) and ( get_next ) ) :
			self.mesg_txt( 'none' )
			return []

		# Convert/standardize the requested time.

		time_req_str = calc_time_str( time )
		time_req_val = calc_time_val( time )
		time_req_epc = calc_time_epc( time )

		# Extract requested date (as a string) and the requested time
		# (as a float indicating seconds since midnight).  Likewise,
		# determine the date of the previous day and the date of the
		# next day.

		date_req_str = time_req_str[0:10]
		scnd_req_val = time_req_val - calc_time_val( date_req_str )

		date_pre_str = ( calc_time_str( time_req_val - 86400. ) )[0:10]
		date_nex_str = ( calc_time_str( time_req_val + 86400. ) )[0:10]

		# Load all the spectra from the requested date.  If the
		# requested time is within "self.buf" seconds of either the
		# previous or next day, load all the spectra from that date as
		# well.

		# Note.  There is no need to check here whether a date has
		#        already been loaded as that's the first thing that
		#        "self.load_date( )" does.

		self.load_date( date_req_str )		

		if ( scnd_req_val <= self.buf ) :
			self.load_date( date_pre_str )

		if ( ( 86400. - scnd_req_val ) <= self.buf ) :
			self.load_date( date_nex_str )

		# If no spectra have been loaded, abort.

		if ( len( self.arr_tag ) == 0 ) :
			self.mesg_txt( 'none' )
			return []

		# Locate the spectrum whose timestamp is closest to the
		# one requested.

		dt  = [ datetime(1970, 1, 1) + timedelta( seconds=tag.epoch ) -
		        self.core.fc_spec['time'] for tag in self.arr_tag      ]

		adt = [ abs( del_t ) for del_t in dt ]

		adt_min = min( adt )

		dt_min = dt[ where( [ del_t == adt_min
		                      for del_t in adt ] )[0][0] ]

		tk = [ a for a in range( len( adt ) ) if adt[a] == adt_min ][0]

#		if ( get_prev ) :
#			tk -= 1
#		if ( get_next ) :
#			tk +=1

		if( ( tk <  0                   ) or
		    ( tk >= len( self.arr_tag ) )    ) :
			self.mesg_txt( 'none' )
			return []

		# Determine how many more PESA-L spectra exist within the next
		# 30 rotations

		num_spec = len( where( [( del_t >= timedelta(seconds=-3.05) and
		                          del_t <= timedelta(seconds=31.*3.05) )
		                         for del_t in dt ] )[0] )

		# If the selected spectrum is not within the the request
		# tolerence, abort.

		if ( ( adt[tk] ).total_seconds() > self.tol ) :# In case of long
		                                               # Data gap  
			self.mesg_txt( 'none' )
			return []

		# Get the PL spectra that lie within this time

		spec = []

		if num_spec == 1 :
			plur = 'spectrum'
		else :
			plur = 'spectra'

		self.mesg_txt( 'load', (str(num_spec) + ' ' + plur + ' found') )

		for n in range( num_spec ) :

			# Extract the spectrum to be returned.

			cdf = self.arr_cdf[self.arr_tag[tk+n].c]
			s   = self.arr_tag[ tk+n ].s

			# Assigning all retrieved data to parameter values

			t_strt   = cdf['sec_beg'][s]

			t_stop   = cdf['sec_end'][s]

	                elev_cen = cdf['the'][s]

			the_del  = cdf['d_the'][s]

	                azim_cen = cdf['phi'][s]

			phi_del  = cdf['d_phi'][s]

	                volt_cen = cdf['nrg'][s]

			volt_del = cdf['d_nrg'][s]

			psd      = cdf['psd'][s]

			spec = spec + [ pl_spec( t_strt=t_strt, t_stop=t_stop,
			                   elev_cen=elev_cen, the_del=the_del,
			                   azim_cen=azim_cen, phi_del=phi_del,
			                   volt_cen=volt_cen, volt_del=volt_del,
			                   psd=psd                           ) ]


		# Request a cleanup of the data loaded into this archive.

		self.cleanup_date( )	

		return spec

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR LOADING ALL SPECTRA FROM DATE-SPECIFIED FILE.
	#-----------------------------------------------------------------------

	def load_date( self, date_str ) :

		# Determine whether or not the requested date has already been
		# loaded.  If it has, abort.

		arr = [ date for date in self.arr_date if date == date_str ]

		if ( len( arr ) > 0 ) :
			return

		# Extract the year, month, and day portions of the "date_str"
		# string.

		str_year = date_str[0:4]
		str_mon  = date_str[5:7]
		str_day  = date_str[8:10]

		# Determine the name of the file that contains data from the
		# requested date.

		fl0 = 'wind-faces_esa_' + \
		       str_year + '-' + str_mon + '-' + str_day + '.idl'

		fl0_path = os.path.join( self.path, fl0 )

		gb = glob( fl0_path )	# returns all files with 
		                        # common names in argument

		# If the file does not exist, say so.

		if ( len( gb ) <= 0 ) :	
			self.mesg_txt( 'fail', date_str )
			return
		else : # Take the last one : gb[-1]
			fl_path = gb[-1]

		# If the file now exists, try to load it; otherwise, abort.

		self.mesg_txt( 'load', date_str )

		if ( os.path.isfile( fl_path ) ) :
			try :
				cdf = readsav( fl_path )
			except :
				self.mesg_txt( 'fail', date_str )
				return
		else :
			self.mesg_txt( 'fail', date_str )
			return

		# Add the CDF object and tags for each spectrum to the arrays.

		c = len( self.arr_cdf )

		self.arr_cdf  = self.arr_cdf  + [ cdf ]	     # arr_cdf and	
		self.arr_date = self.arr_date + [ date_str ] # arr_date of
		                                             # same size

		n_spec = len( cdf['sec_beg'] )
		self.arr_tag = self.arr_tag + [ pl_tag( c=c, s=s,
		                                epoch=cdf['sec_beg'][s] )
		                                for s in range( n_spec ) ]

		self.arr_tag = sorted( self.arr_tag, key=attrgetter('epoch') )




		# Request a clean-up of the files in the data directory.

		#self.cleanup_file( )



	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR CLEANING UP THIS ARCHIVE.
	#-----------------------------------------------------------------------

	def cleanup_date( self ) :

		if ( len(self.arr_date) <= 0 ) :
			return

		# How to get the entire list of arr_cdf, arr_data, arr_tag for 
		#all downloded data. ...Showing for only the requested one ...

		self.arr_tag = [ tag for tag in self.arr_tag if tag.c != 0 ]

		self.arr_cdf  = self.arr_cdf[1:]
		self.arr_date = self.arr_date[1:]

		for t in range( len( self.arr_tag ) ) :
			self.arr_tag[t].c -= 1

		# In case there are multiple dates to be removed; \
		#special case time 12 a.m.
		
		self.cleanup_date()

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR SENDING INFORMATIONAL MESSAGES TO THE USER.
	#-----------------------------------------------------------------------

	def mesg_txt( self, mesg_typ, mesg_obj='' ) :

		# If this object is not associated with an instance of the
		# Janus "core", or if messaging has not been requested, abort.

		if ( ( self.core is None ) or ( not self.verbose ) ) :

			return

		# Emit a message signal (on behalf of the core) containing the
		# message parameters.

		self.core.emit( SIGNAL('janus_mesg'),
		                'pl', mesg_typ, mesg_obj )
