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
from PyQt4.QtGui import QTextCursor, QScrollBar

# Load the modules for displaying text output.

from janus_format_TextEdit import format_TextEdit
from PyQt4.QtGui import QTextEdit

from janus_widget_mom_res import widget_mom_res


################################################################################
## DEFINE THE "widget_mom_res" CLASS TO CUSTOMIZE "format_TextEdit".
################################################################################

class widget_mom_res_table( format_TextEdit ) :

	#-----------------------------------------------------------------------
	# DEFINE THE INITIALIZATION FUNCTION.
	#-----------------------------------------------------------------------

	def __init__( self, core ) :

		# Inherit all attribues of an instance of "format_TextEdit".

		super( widget_mom_res_table, self ).__init__( )

		# Store the Janus core.

		self.core = core

		# Prepare to respond to signals received from the Janus core.

		self.connect( self.core, SIGNAL('janus_rset'), self.resp_rset )
		self.connect( self.core, SIGNAL('janus_chng_opt'),
		                                            self.resp_chng_opt )
		self.connect( self.core, SIGNAL('janus_chng_mfi'),
		                                            self.resp_chng_mfi )
		self.connect( self.core, SIGNAL('janus_chng_mom_res'),
		                                        self.resp_chng_mom_res )

		# Set this text editor as read only (for the user).

		self.setReadOnly( True )
		self.setLineWrapMode( 0 )

		# Populate this text editor with general information about this
		# spectrum's MFI data.

		self.make_txt( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR GENERATING THE TEXT FOR THIS TEXT AREA.
	#-----------------------------------------------------------------------

	def make_txt( self ) :

		# Clear the text area (just in case there's some text already
		# there).

		self.clear( )

		# If the moments analysis has failed or has not been performed,
		# return.

		if ( self.core.mom_res is None ) :
			return

		self.fc_text = None
		self.fc_text = widget_mom_res( self.core )
		self.fc_string = self.fc_text.toPlainText( )
		print self.fc_string

		# Generate the table for displaying FC and PL data

		self.insertHtml('<table style="width:100%"><tr><th>FC</th><th>PL</th></tr><tr><td>{}</td><td>{}</td></tr></table>'.format(self.fc_string, 'data') )

		# Scroll to the top of the text area.

		self.moveCursor( QTextCursor.Start )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "rset" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_rset( self ) :

		# Reset the text area.

		self.clear( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "chng_opt" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_chng_opt( self ) :

		# Regenerate the text in the text area.

		self.make_txt( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "chng_mfi" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_chng_mfi( self ) :

		# Replace the text in the text area.

		self.make_txt( )

	#-----------------------------------------------------------------------
	# DEFINE THE FUNCTION FOR RESPONDING TO THE "chng_mom_res" SIGNAL.
	#-----------------------------------------------------------------------

	def resp_chng_mom_res( self ) :

		# Replace the text in the text area.

		self.make_txt( )
