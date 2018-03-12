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

from PyQt4.QtCore import Qt, QPointF, SIGNAL
from PyQt4.QtGui import QGridLayout, QWidget, QLabel

# Load the modules necessary for plotting.

from pyqtgraph import AxisItem, GraphicsLayoutWidget, LabelItem, mkBrush, \
                      mkPen, PlotDataItem, TextItem

from janus_event_ViewBox import event_ViewBox

# Load the module necessary handling step functions.

from janus_step import step

# Load the necessary "numpy" array modules and numeric-function modules.

from numpy import amax, amin, array, ceil, floor, log10, sqrt, tile, where

# Load the necessary threading modules.

from threading import Thread
from janus_thread import n_thread, thread_chng_mom_sel, thread_chng_nln_sel

# Load the modules necessary handling dates and times.

from datetime import datetime, timedelta

# Load the module for TESTING joint

from scipy.io import readsav
from janus_pl_spec import pl_spec
from numpy import zeros

################################################################################
## DEFINE THE "widget_fc_cup" CLASS TO CUSTOMIZE "QWidget" FOR Wind/FC PLOTS.
################################################################################

class widget_pl_grid( QWidget ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core, n,
	              n_plt_x=None, n_plt_y=None, n_plt=None ) :

		# Inherit all attributes of an instance of "QWidget".

		super( widget_pl_grid, self ).__init__( )

		# Initialize the counter of repaint events for this widget as
		# well as a maximum value for this counter.

		# Note.  For some reason, adjusting the individual plots to have
		#        uniform sizes is difficult to achieve before the widget
		#        is rendered.  Thus, once a paint event occurs, the
		#        "self.paintEvent( )" function picks it up and makes a
		#        call to "self.ajst_grd( )".  This counter and its
		#        maximum value are used ensure that "self.paintEvent( )"
		#        makes such a call only in response to the intial few
		#        painting (so as to prevent an infinite loop).

		# Note.  The first paint seems to be a "dummy" of some sort.
		#        Whatever the case, "self.n_paint_max = 1" seems to
		#        generally be insufficient.

		self.n_painted     = 0
		self.n_painted_max = 3

		self.core = core
		self.n = n
		self.t = []
		self.delta_t = []
		self.time_label = QLabel( )
		self.time_label.setAlignment( Qt.AlignCenter )

		# Prepare to respond to signals received from the Janus core.

		self.connect( self.core, SIGNAL('janus_rset'), self.resp_rset )
		self.connect( self.core, SIGNAL('janus_chng_pl_spc'),
		                                            self.resp_chng_pl_spc )
		#TODO add more signals

		# Assign (if not done so already) and store the shape of the
		# plot-grid array.

		self.n_plt_x = 5 if ( n_plt_x is None ) else n_plt_x
		self.n_plt_y = 5 if ( n_plt_y is None ) else n_plt_y

		if ( n_plt is None ) :
			self.n_plt = self.n_plt_x * self.n_plt_y

		# Initizalize the pens, brushes, and fonts used by this widget.

		self.pen_plt   = mkPen( color='k' )
		self.pen_hst   = mkPen( color='k' )
		self.pen_pnt_c = mkPen( color='k' )
		self.pen_pnt_y = mkPen( color='k' )
		self.pen_pnt_r = mkPen( color='k' )
		self.pen_crv_b = mkPen( color='b' )
		self.pen_crv_g = mkPen( color='g' )

		self.bsh_pnt_c = mkBrush( color='c' )
		self.bsh_pnt_y = mkBrush( color='y' )
		self.bsh_pnt_r = mkBrush( color='r' )

		self.fnt = self.core.app.font( )

		# Set the maximum number of velocity channels and the maximum
		# number of ion species.

		self.n_k   = 14
		self.n_ion = self.core.nln_n_pop

		# Initialize the widget and it's plot's.

		self.init_plt( )

		# Populate the plots with the histograms (and labels), the
		# selection points, and the fit curves.

		self.make_hst( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR INITIALIZING THE WIDGET AND ITS PLOTS.
	#-----------------------------------------------------------------------

	def init_plt( self ) :

		# Initialize the "GraphicsLayoutWidget" for this widget.  This
		# will allow a grid of "GraphicsItem" objects, which will
		# include the plots themselves, the axes, and the axis labels.

		# Note.  The "QGridLayout" object given to this widget as its
		#        layout is essentially a dummy.  I tried to just having
		#        this widget be an extention of "GraphicsLayoutWidget"
		#        (i.e., having it inheret that type), but I couldn't get
		#        it to display anything at all.

		self.setLayout( QGridLayout( ) )

		self.grd = GraphicsLayoutWidget( )
		self.grd.setBackground( 'w' )
		self.layout( ).addWidget( self.time_label )
		self.layout( ).addWidget( self.grd )

		self.layout().setContentsMargins( 0, 0, 0, 0 )

		# Initialize the text for the x- and y-axis labels.  Then,
		# create the labels themselves and add them to the grid.

		self.txt_axs_x = 'Projected Proton Inflow Velocity [km/s]'
		self.txt_axs_y = u'Phase-space Density [cm\u00AF\u00B3/(km/s)\u00B3]'

		if ( self.core.app.res_lo ) :
			size =  '8pt'
		else :
			size = '10pt'

		self.lab_axs_x = LabelItem( self.txt_axs_x, angle=0  ,
		                            color='b', size=size       )
		self.lab_axs_y = LabelItem( self.txt_axs_y, angle=270, 
		                            color='b', size=size       )

		self.grd.addItem( self.lab_axs_x, self.n_plt_y + 1, 2,
		                                  1, self.n_plt_x      )
		self.grd.addItem( self.lab_axs_y, 0, 0,
		                                  self.n_plt_y, 1      )

		# Initialize the arrays that will contain the individual axes,
		# plots, and plot elements (i.e., the histograms, fit curves,
		# labels, and selection points).

		self.plt = tile( None, [ self.n_plt_y, self.n_plt_x ] )

		self.axs_x = tile( None, self.n_plt_x )
		self.axs_y = tile( None, self.n_plt_y )

		self.hst = tile( None, [ self.n_plt_y, self.n_plt_x ] )
		self.lbl = tile( None, [ self.n_plt_y, self.n_plt_x ] )

		# Initialize the scale-type for each axis, then generate the
		# (default) axis-limits and adjusted axis-limits.

		self.log_x = False
		self.log_y = True

		self.make_lim( )

		# Create, store, and add to the grid the individual axes: first
		# the horizontal and then the vertical.

		for i in range( self.n_plt_x ) :

			self.axs_x[i] = AxisItem( 'bottom', maxTickLength=5 )
			self.axs_x[i].setLogMode( self.log_x )
			self.axs_x[i].setRange( self.x_lim[0], self.x_lim[1] )
			self.axs_x[i].setTickFont( self.fnt )

			if ( self.core.app.res_lo ) :
				self.axs_x[i].setHeight( 10 )
			else :
				self.axs_x[i].setHeight( 20 )

			self.grd.addItem( self.axs_x[i], self.n_plt_y, i + 2 )

		for j in range( self.n_plt_y ) :

			self.axs_y[j] = AxisItem( 'left', maxTickLength=5 )
			self.axs_y[j].setLogMode( self.log_y )
			self.axs_y[j].setRange( self.y_lim[0], self.y_lim[1] )
			self.axs_y[j].setTickFont( self.fnt )

			if ( self.core.app.res_lo ) :
				self.axs_y[j].setWidth( 32 )
			else :
				self.axs_y[j].setWidth( 40 )

			self.grd.addItem( self.axs_y[j], j, 1 )

		# Create, store, and add to the grid the individual plots.
		# Likewise, create, store, and add to each plot a label.

		for t in range( self.n_plt_y ) :

			for p in range( self.n_plt_x ) :

				# Compute the plot number of this plot.

				d = p + ( t * self.n_plt_x )

				# If creating this plot would exceed the
				# specified number of plots, don't create it.

				if ( d >= self.n_plt ) :
					continue

				# Create and store this plot, adjust its limits,
				# and add it to the grid.

				self.plt[t,p] = event_ViewBox( self,
				                          border=self.pen_plt,
				                          enableMouse=False,
				                          enableMenu=False     )

				self.plt[t,p].setRange( xRange=self.x_lim,
				                        yRange=self.y_lim,
				                        padding=0.         )

				self.grd.addItem( self.plt[t,p], t, p + 2 )

				# Create and store an (empty) label and add it
				# to this plot.

				self.lbl[t,p] = TextItem( anchor=(1,0) )

				self.lbl[t,p].setFont( self.fnt )

				self.plt[t,p].addItem( self.lbl[t,p] )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR GENERATING AXIS-LIMITS (AND ADJUSTED LIMITS).
	#-----------------------------------------------------------------------

	def make_lim( self ) :

		# If no spectrum has been loaded, use the default limits;
		# otherwise, use the spectral data to compute axis limits.

		if ( self.core.pl_spec_arr is None ) :

			self.domain = [ 300. , 900. ]
			self.range = [ 1.e-10,  1.e-5 ]

		else :

			self.domain = [self.core.pl_spec_arr[self.n]['vel_strt'][0],
			               self.core.pl_spec_arr[self.n]['vel_stop'][-1]]

			arr_psd_flat = [self.core.pl_spec_arr[self.n]['psd_flat'][i] for i
			          in (where(array(self.core.pl_spec_arr[self.n]['psd_flat'])
			                                             != 0.)[0])]

			self.range = [ min(arr_psd_flat), max(arr_psd_flat)  ]

			# Note: psd values are less than 1

			if ( self.log_y ) :
				self.range[1] = self.range[1] ** 0.9
			else :
				self.range[1] += 0.1 * ( self.range[1] -
				                         self.range[0]   )

		# Compute the "adjusted limits" for each axis.

		if ( self.log_x ) :
			self.x_lim = [ log10( x ) for x in self.domain ]
		else :
			self.x_lim = self.domain

		if ( self.log_y ) :
			self.y_lim = [ log10( y ) for y in self.range ]
		else :
			self.y_lim = self.range

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR CREATING THE PLOTS' HISTOGRAMS (AND LABELS).
	#-----------------------------------------------------------------------

	def make_hst( self ) :

		# If no spectrum has been loaded, clear any existing histograms
		# and abort.

		if ( self.core.pl_spec_arr is None ) : return

		# Generate the timestamp label

		self.t = self.t + [self.core.pl_spec_arr[self.n]['time'][0]]

		self.t_0 = self.core.fc_spec['time']

		self.delta_t = self.delta_t + [(self.t[-1]-self.t_0).total_seconds( )]

		self.time_label.setText( str(self.t[-1]) + '        ' + u'\u0394t = {}'.format(
		                                            self.delta_t[-1] ) )

		# Use the spectral data to compute new axis-limits.

		self.make_lim( )

		for p in range( self.n_plt_x ) :
			 self.axs_x[p].setRange( self.x_lim[0], self.x_lim[1] )

		for t in range( self.n_plt_y ) :
			 self.axs_y[t].setRange( self.y_lim[0], self.y_lim[1] )

		# Histograms are broken down by phi horizontally and
		# theta vertically

		for p in range( self.core.pl_spec_arr[self.n]['n_phi'] ):

			for t in range ( self.core.pl_spec_arr[self.n]['n_the'] ):

				# If this plot does not exist, move onto
				# the next one.

				if ( self.plt[t,p] is None ) :
					continue

				#-----------------------------#
				#---DATA GENERATION SECTION---#
				#-----------------------------#

				# Generate a step function for the
				# look direction associated with this widget.

				self.stp = array( [step( self.core.pl_spec_arr[self.n]['vel_cen'],
				                        self.core.pl_spec_arr[self.n]['vel_del'],
				                        self.core.pl_spec_arr[self.n]['psd'][t][p])])

				# Calculate the points to be plotted from the
				# step function

				stp_pnt = array( [ array( datum.calc_pnt( 
				                     lev_min=self.range[0]/2.) )
				                   for datum in self.stp     ] )

				self.x_set = stp_pnt[:,0][0]
				self.y_set = stp_pnt[:,1][0]

				# If plotting log(y) and there are any psd
				# values of zero, replace those points with an
				# arbitrary minimum y value

				if ( self.log_y ) :
					y_min = self.range[0]/2.
					self.y_lim[0] = log10(y_min)
					self.y_set = [ max( y, y_min ) for y
					                         in self.y_set ]

				# If generating a log plot, take the log of the
				# points to be plotted

				self.x_pnts = log10( self.x_set ) if ( self.log_x ) else \
				                     self.x_set
				self.y_pnts = log10( self.y_set ) if ( self.log_y ) else \
				                     self.y_set

				#---------------------------------#
				#---GRAPHICS GENERATION SECTION---#
				#---------------------------------#

				# If a histogram already exists for this plot,
				# remove and delete it.

				if ( self.hst[t,p] is not None ) :
					self.plt[t,p].removeItem(self.hst[t,p])
					self.hst[t,p] = None

				# Clear this plot's label of text.

				self.lbl[t,p].setText( '' )

				# Adjust this plot's limits and then move it's
				# label in response.

				self.plt[t,p].setRange( xRange=self.x_lim,
				                        yRange=self.y_lim,
				                        padding=0.         )

				self.lbl[t,p].setPos( self.x_lim[1],
				                      self.y_lim[1]  )

				# Update this plot's label with appropriate text
				# indicating the pointing direction.


				elev = round( self.core.pl_spec_arr[self.n]['elev_cen'][t] )
				azim = round( self.core.pl_spec_arr[self.n]['azim_cen'][p] )

				txt = ( u'({0:+.0f}\N{DEGREE SIGN}, ' + 
				        u'{1:+.0f}\N{DEGREE SIGN})'     ).format(
				                                          elev, azim )

				self.lbl[t,p].setText( txt, color=(0,0,0) )

				# Generate the histogram for the data from this
				# look direction and display it in the plot.

				self.hst[t,p] = PlotDataItem( self.x_pnts,
				                              self.y_pnts,
				                              pen=self.pen_hst )

				self.plt[t,p].addItem( self.hst[t,p] )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESETTING THE PLOTS' HISTOGRAMS (AND LABELS).
	#-----------------------------------------------------------------------

	def rset_hst( self, rset_lbl=False ) :

		self.time_label.clear( )
		self.t = []
		self.delta_t = []

		# For each plot that exists in the grid, remove and delete it's
		# histogram.  Likewise, if requested, empty it's label (but
		# still leave the label itself intact).

		for t in range( self.n_plt_y ) :

			for p in range( self.n_plt_x ) :

				# If the plot does not exist, move onto the the
				# next one.

				if ( self.plt[t,p] is None ) :
					continue

				# If a histogram exists for this plot, remove
				# and delete it.

				if ( self.hst[t,p] is not None ) :
					self.plt[t,p].removeItem(
					                         self.hst[t,p] )
					self.hst[t,p] = None

				# If requested, reset this plot's label text to
				# the empty string.

				if ( rset_lbl ) :
					self.lbl[t,p].setText( '',
					                       color=(0,0,0) )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "rset" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_rset( self ) :

		# Clear the plots of all their elements.

		self.rset_hst( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "chng_pl_spc" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_chng_pl_spc( self ) :

		# Clear the plots of all their elements and regenerate them.

		self.rset_hst( )

		#self.make_hst( )
