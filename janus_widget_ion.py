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

from PyQt4.QtGui import QTabWidget

from janus_widget_ion_fc import widget_fc

from janus_widget_ion_pl import widget_pl

################################################################################
## DEFINE THE "widget_ion" CLASS TO CUSTOMIZE "QTabWidget" FOR Wind/FC AND
## Wind/PESA-L PLOTS.
################################################################################

class widget_ion( QTabWidget ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core ) :

		# Inherit all attributes of an instance of "QTabWidget".

		super( widget_ion, self ).__init__( )

		# Store the Janus core.

		self.core = core

		# Create one instance of "widget_fc" and one instance of
		# "widget_pl and add each as a tab.

		self.wdg_fc = widget_fc( core=self.core )
		self.wdg_pl = widget_pl( core=self.core, n=1 )

		self.addTab( self.wdg_fc, 'Faraday Cup' )
		self.addTab( self.wdg_pl,  'PESA-Low' )

		self.connect( self.core, SIGNAL('janus_pl_reset'), self.pl_reset )

	def pl_reset( self ) :
		self.removeTab(1)
		self.wdg_pl = widget_pl( core=self.core, n=1 )
		self.addTab( self.wdg_pl,  'PESA-Low' )
		

