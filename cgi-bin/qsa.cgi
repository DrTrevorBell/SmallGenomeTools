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

form = cgi.FieldStorage()

version = "Quality Score Analyzer"

def closePage():
	print '<hr>'
	print '<tt><p align=right>Version 0.1 (July 2010)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'Quality Score Analyzer\t' )
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

inFilename = form["infile"].filename

if len(inFilename) == 0:
	badFile('Specify a chromatogram file!', "")

inFile = form["infile"].file

print '<h1> <style="vertical-align:middle">Quality Score Analyzer Results</h1><br>'

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'QSA' + randomToken + '/'
os.makedirs(climbMini.TEMPFOLDER)

tempF = climbMini.TEMPFOLDER + inFilename
temp = open(tempF, 'wb')
temp.write(inFile.read())
temp.close()
inFile.close()

X = 5
Y = 20

import ABIFReader
reader = ABIFReader.ABIFReader(tempF)
qualTemp = reader.getData('PCON')	# returns a string
base = reader.getData('PBAS')		# returns a string
seqID = reader.getData('SMPL')		# sample name (FASTA identifier)
reader.close()

qual = []
for i in qualTemp:
	qual.append(ord(i))		# store quality scores as integers in a list

htmlPrint('<b>Input file:</b> %s' % inFilename)
htmlPrint('<b>Sequence ID:</b> %s' % seqID)
htmlPrint('<b>Sequence Length:</b> %i' % len(base))
htmlPrint('<b>Minimum Quality Score:</b> %i' % min(qual))
htmlPrint('<b>Maximum Quality Score:</b> %i' % max(qual))

import rpy

boxFile = climbMini.TEMPFOLDER + 'qsBox.png'
rpy.r.png(boxFile)
rpy.r.boxplot(qual, ylim=(0,max(qual)+5), las=3, ylab="Quality Score")
rpy.r('dev.off()')

densityFile = climbMini.TEMPFOLDER + 'qsDensity.png'
rpy.r.png(densityFile)
rpy.r.plot(rpy.r.density(qual), xlab="Quality Score", ylab="Density", xlim=(-5, max(qual)+10), type="l")
rpy.r('dev.off()')

print '<img src=%s><img src=%s><br>' % (boxFile[8:], densityFile[8:])

htmlPrint("<h2>Quality Score Heat Map</h2>")

print '<style type="text/css" rel="stylesheet">'
print '<!--'
print '.nowrap {font-family: monospace; font-size: 12px; font-weight: normal; background: #000000; color: #FFFFFF; white-space: nowrap; border:2px solid black; overflow-y:hidden; overflow-x:scroll; width:760px}'
print '-->'
print '</style>'

print '<div class="nowrap">'

basesPerLine = 12

RED =     '<font color="red">'
GREEN =   '<font color="green">'
YELLOW =  '<font color="yellow">'
BLUE =    '<font color="blue">'
MAGENTA = '<font color="magenta">'
CYAN =    '<font color="cyan">'
WHITE =   '<font color="white">'
RESET =   '</font>'

RED_ =     '<font style="background: red" color="black">'
GREEN_ =   '<font style="background: green" color="black">'
YELLOW_ =  '<font style="background: yellow" color="black">'
BLUE_ =    '<font style="background: blue" color="black">'
MAGENTA_ = '<font style="background: magenta" color="black">'
CYAN_ =    '<font style="background: cyan" color="black">'
WHITE_ =   '<font style="background: white" color="black">'

out = ''
for i in range(len(qual)):
	b = base[i]
	q = qual[i]
	if b.upper() in ['A','C','G','T']:
		if q < 10:
			out += RED
		elif q < 20:
			out += YELLOW
		elif q < 30:
			out += GREEN
		elif q < 40:
			out += BLUE
		elif q < 50:
			out += MAGENTA
		elif q < 60:
			out += CYAN
		else:
			out += WHITE
	else:
		if q <= 10:
			out += RED_
		elif q <= 20:
			out += YELLOW_
		elif q <= 30:
			out += GREEN_
		elif q <= 40:
			out += BLUE_
		elif q <= 50:
			out += MAGENTA_
		elif q <= 60:
			out += CYAN_
		else:
			out += WHITE_


	out += '%04i:%02i%s' % (i+1, q, b) + RESET + ' '		# prevent underlined / inverse spaces
	if (i+1) % basesPerLine == 0:	# i starts at 0; 0 % basesPerLine  = 0
		htmlPrint(out)
		out = ''
		partial = False
	else:
		partial = True

if partial:
	htmlPrint(out)

print '</p></div><br><br>'

closePage()

