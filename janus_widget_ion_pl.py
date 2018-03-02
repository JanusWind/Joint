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

		# Create two instances of "widget_fc_cup" (one for each Faraday
		# cup) and add each as a tab.

		self.wdg_pl_grid1 = widget_pl_grid( core=self.core,
		                              n_plt_x=5, n_plt_y=5 )
#		self.wdg_pl_cont1 = widget_pl_cont( core=self.core,
#		                              n_plt_x=n_plt_x, n_plt_y=n_plt_y )

		self.addTab( self.wdg_pl_grid1, 'PESA-Low Grid 1' )
#		self.addTab( self.wdg_pl_cont1, 'PESA-Low Countour 1'








