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

# Load the modules necessary for the graphical interface.

from PyQt4.QtGui import QTabWidget, QWidget

from janus_widget_ion_pl_grid import widget_pl_grid

#from janus_widget_ion_pl_cont import widget_pl_cont

################################################################################
## DEFINE THE "widget_fcspec" CLASS TO CUSTOMIZE "QTabWidget" FOR Wind/FC PLOTS.
################################################################################

class widget_pl( QTabWidget ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core ) :

		# Inherit all attributes of an instance of "QTabWidget".

		super( widget_pl, self ).__init__( )

		# Store the Janus core.

		self.core = core

		# Prepare to respond to signals received from the Janus core.

		self.connect( self.core, SIGNAL('janus_reset'), self.resp_rset )
		self.connect( self.core, SIGNAL('janus_chng_pl_spc'),
		                                         self.resp_chng_pl_spc )

		# Initialize the widget array

		self.wdg_arr = []

		for n in range(8) :

			wdg = widget_pl_grid( self.core, n,
			                          n_plt_x=5, n_plt_y=5 )

			self.wdg_arr = self.wdg_arr + [ wdg ]

#			self.wdg_pl_cont = widget_pl_cont( core=self.core,
#		                              n_plt_x=n_plt_x, n_plt_y=n_plt_y )

			self.addTab( wdg, 'P-L Grid {}'.format(n+1) )

		self.hide_tabs( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR REMOVING THE TABS.
	#-----------------------------------------------------------------------

	def hide_tabs( self ) :

		return

		# Hide all instances of widget_arr tabs in the GUI

		for wdg in self.wdg_arr :

			wdg.hide( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "rset" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_rset( self ) :

		# Clear any current tab widgets

		self.hide_tabs( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "chng_pl_spec" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_chng_pl_spc( self ) :

		# Clear any current tab widgets

		self.hide_tabs( )

		# Create each instance of "widget_pl_grid" and each
		# instance of "widget_pl_cont" and add each as a tab.

		for n in range( len( self.core.pl_spec_arr ) ) :

			self.wdg_arr[n].show( )
