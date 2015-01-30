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

#!/usr/bin/python -u

# Output buffering disabled by "-u" above

import datetime

starttime = datetime.datetime.now()

import cgi, sys, socket, re, os, time

import climbMini

import cgitb
cgitb.enable()

form = cgi.FieldStorage()

version = "Automatic Contig Generator Tool (ACGT)"

def closePage():
	print '<hr>'
	print '<tt><p align=right>Version 0.2 (May 2010)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'ACGT\t' )
	f.write( 'Random Token %s\t' % (randomToken) )
	f.write( 'Served by %s at %s\t' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0]) )
	f.write( 'Run completed %s\t' % (endtime.strftime("%Y-%m-%d %H:%M")) )
	f.write( '%03.6f seconds\t' % ((endtime - starttime).microseconds / float(1000000)) )
	f.write( '\n')
	f.close()

def badFile(errorText, explainText):
	print "<h1>%s</h1>" % (version)
	print "<h2>%s</h2>" % (errorText)
	print explainText + '<br>'
	closePage()
	sys.exit()

def htmlPrint(pp, lines=1):
	print pp
	for i in range(lines):
		print "<br>"
# ---------------

print 'Content-Type: text/html'
print

print '<!DOCTYPE html'
print '	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
print '	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">'
print '<head><title>HVDRP %s</title>' % (version)
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head><body>'

randomToken = climbMini.randomStamp()

symbolEqual = form["symbolEqual"].value
symbolMismatch = form["symbolMismatch"].value
symbolGap = form["symbolGap"].value
fontSize = form["fontSize"].value

htmlPrint("<h1>Contig Generated</h1>")

forwardFile = form["forwardfile"].file
reverseFile = form["reversefile"].file
forwardFilename = form["forwardfile"].filename
reverseFilename = form["reversefile"].filename

if len(forwardFilename) == 0 or len(reverseFilename) == 0:
	badFile('Error: Missing File!', '')

climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'ACGT-' + randomToken + '/'

os.makedirs(climbMini.TEMPFOLDER)

tempF = climbMini.TEMPFOLDER + forwardFilename
temp = open(tempF, 'wb')
temp.write(forwardFile.read())
temp.close()
forwardFile.close()

tempR = climbMini.TEMPFOLDER + reverseFilename
temp = open(tempR, 'wb')
temp.write(reverseFile.read())
temp.close()
forwardFile.close()

S1 = climbMini.Sequence()
S2 = climbMini.Sequence()

S1.load(tempF, ABIF=True)
S2.load(tempR, ABIF=True)

X = climbMini.contig(S1, S2)

htmlPrint('<b>Forward chromatogram:</b> ' + forwardFilename)
htmlPrint('<b>Reverse chromatogram:</b> ' + reverseFilename, lines=2)
htmlPrint('<i>Right-click the <img src="/res/download.png"> image to download a FASTA file of the data in each box. The characters in white above the alignment indicate the level of consensus: a </i><b>%s</b><i> character indicates that the bases in the alignment are equal; a </i><b>%s</b><i> character indicates that the bases in the alignment are mismatched; a </i><b>%s</b><i> character indicates that one of the bases in the alignment is a gap.</i>' % (symbolEqual, symbolMismatch, symbolGap), lines=2)

print '<style type="text/css" rel="stylesheet">'
print '<!--'
print '.nowrap {font-family: monospace; font-size: %spx; font-weight: normal; background: #000000; white-space: nowrap; border:1px solid black; overflow-y:hidden; overflow-x:scroll;}' % (fontSize)
print '-->'
print '</style>'

newline = '\n'

contigFile = climbMini.TEMPFOLDER + X[0][0][1:] + '.fasta'
f = open(contigFile, 'w')
f.write(X[0][0] + newline)
f.write(X[0][1])
f.close()

alignFile = climbMini.TEMPFOLDER + X[2][0][1:] + '.fasta'
f = open(alignFile, 'w')
f.write(X[2][0] + newline)
f.write(X[2][1] + newline)
f.write(X[3][0] + newline)
f.write(X[3][1])
f.close()

htmlPrint('<a href="%s"><img src="/res/download.png" style="border-style: none"></a><b><font color="#BBBB00">Alignment</font></b> &#8226; <a href="%s"><img src="/res/download.png" style="border-style: none"></a><b><font color="#FF00FF">Contig</font></b> ' % (alignFile[8:], contigFile[8:]))

print '<div class="nowrap">'

conLen = len(X[0][1])
vis = ''
mismatches = 0
for i in range(conLen):
	if X[2][1][i] == '-' or X[3][1][i] == '-':
		vis += symbolGap
	elif X[2][1][i] == X[3][1][i]:
		vis += symbolEqual
	else:
		vis += symbolMismatch
		mismatches += 1

# LEVEL OF CONSENSUS
print '<font color="#FFFFFF">%s</font><br>' % (vis)

# ALIGNMENT
print '<font color="#FFFF00">'
htmlPrint(X[2][1])
htmlPrint(X[3][1])
print '</font>'

# CONTIG
print '<font color="#FF00FF">'
htmlPrint(X[0][1])
print "</font>"

print '</p></div>'
htmlPrint('')

gaps = X[2][1].count('-') + X[3][1].count('-')
equal = conLen - mismatches - gaps

htmlPrint('<b>Contig Quality Score:</b> %5.4f' % (X[1]), lines=2)

print '<table border="1">'
print '<tr><td><b>Length</b></td><td><b>Equal (%s)</b></td><td><b>Mismatch (%s)</b></td><td><b>Gap (%s)</b></td></b></tr>' % (symbolEqual, symbolMismatch, symbolGap)
print '<tr><td>%i</td><td>%i</td><td>%i</td><td>%i</td></tr>' % (conLen, equal, mismatches, gaps)
print '<tr><td>%3.2f%%</td><td>%3.2f%%</td><td>%3.2f%%</td><td>%3.2f%%</td></tr>' % (conLen/float(conLen)*100, equal/float(conLen)*100, mismatches/float(conLen)*100, gaps/float(conLen)*100)
print '</table><br>'

inputFile = climbMini.TEMPFOLDER + X[4][0][1:] + '.fasta'
f = open(inputFile, 'w')
f.write(X[4][0] + newline)
f.write(X[4][1] + newline)
f.write(X[5][0] + newline)
f.write(X[5][1])
f.close()
htmlPrint('<a href="%s"><img src="/res/download.png" style="border-style: none"></a><b><font color="#00BBBB">Trimmed Input Sequences</font></b>' % (inputFile[8:]))
print '<div class="nowrap">'
print '<p style="color:#00FFFF">'
htmlPrint(X[4][0])
htmlPrint(X[4][1])
htmlPrint(X[5][0])
htmlPrint(X[5][1])
print '</p></div>'
htmlPrint('')

closePage()

