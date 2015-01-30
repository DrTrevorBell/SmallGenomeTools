# Small Genome Tools
# Copyright (C) 2015 University of the Witwatersrand, Johannesburg, South Africa
# Author: Dr Trevor G. Bell, TrevorGrahamBell@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#!/usr/bin/python

import datetime

starttime = datetime.datetime.now()

import cgi, os, sys, socket

import cgitb
cgitb.enable()

import climbMini

form = cgi.FieldStorage()

version = "<h1>Automatic Alignment Clean-up Tool (AACT)</h1>"

def closePage():
	print '<hr>'
	print '<tt><p align=right>Version 0.1 Beta (July 2010)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'AACT\t' )
	f.write( 'Served by %s at %s\t' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0]) )
	f.write( 'Run completed %s\t' % (endtime.strftime("%Y-%m-%d %H:%M")) )
	f.write( '%03.6f seconds\t' % ((endtime - starttime).microseconds / float(1000000)) )
	f.write( '\n')
	f.close()

def badFile(errorText, explainText):
	print version
	print "<h2>%s</h2>" % (errorText)
	print explainText
	closePage()
	# infile.close()
	sys.exit()

print 'Content-Type: text/html'
print

print '<!DOCTYPE html'
print '	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
print '	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">'
print '<head><title>HVDRP Sequence Merger</title>'
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head><body>'

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file!', '')

eliminate = form.has_key("eliminate")
disambiguate = form.has_key("disambiguate")
bluntLeft = form.has_key("bluntleft")
bluntRight = form.has_key("bluntright")

gapThresholdHere = int(form["gapthresholdpercent"].value.rstrip()) / float(100)

if (eliminate) and (gapThresholdHere == 0):
	badFile('Specify gap threshold percent')

if not eliminate and not disambiguate:
	badFile('No action specified!', '')

print '<h1> <style="vertical-align:middle">Automatic Alignment Clean-up Tool (AACT) Results</h1><br>'

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'AACT-' + randomToken + '/'
os.makedirs(climbMini.TEMPFOLDER)

tempF = climbMini.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

S = climbMini.Sequence()
S.load(tempF)

if bluntLeft or bluntRight:
	bluntResult = S.blunt(left=bluntLeft, right=bluntRight)
	if bluntLeft:
		print 'Columns removed from start of sequence (blunt): %i<br>' % bluntResult[0]
	if bluntRight:
		print 'Columns removed from end of seuqence (blunt): %i<br>' % bluntResult[1]
	print '<br>'

if eliminate:
	result = S.eliminateGapColumns(gapThreshold=gapThresholdHere, verbose=True)

	for line in result:
		print '%s<br>' % line

	print '<br>'

if disambiguate:
	result = S.disambiguateColumns(verbose=True)

	for line in result:
		print '%s<br>' % line

	print '<br>'

S.save(S.fileName, overWrite=True)

print('<a href="%s"><img src="/res/download.png" style="border-style: none"></a> Right-click icon and select "Save As" to download file in FASTA format<br><br>' % (S.fileName[8:]))

S.unload()

closePage()

