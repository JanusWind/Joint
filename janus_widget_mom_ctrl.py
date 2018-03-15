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

# Load the modules necessary for the graphical interface.

from PyQt4.QtGui import QTabWidget

from janus_widget_mom_ctrl_fc import widget_mom_ctrl_fc
from janus_widget_mom_ctrl_pl import widget_mom_ctrl_pl

################################################################################
## DEFINE THE "widget_mom" CLASS TO CUSTOMIZE "QTabWidget" FOR MOMENTS ANALYSIS.
################################################################################

class widget_mom_ctrl( QTabWidget ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core ) :

		# Inherit all attributes of an instance of "QTabWidget".

		super( widget_mom_ctrl, self ).__init__( )

		# Store the Janus core.

		self.core = core

		# Intialize this widget's sub-widgets and add them as tabs.

		self.wdg_ctrl_fc = widget_mom_ctrl_fc( self.core )
		self.wdg_ctrl_pl = widget_mom_ctrl_pl(  self.core )

		self.addTab( self.wdg_ctrl_fc, 'FC' )
		self.addTab( self.wdg_ctrl_pl, 'PL' )
