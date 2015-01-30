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

version = "<h1>""Purdy"" HBV Serotyper</h1>"

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
	f.write( 'Serotyper\t' )
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
	sys.exit()

print 'Content-Type: text/html'
print

print '<!DOCTYPE html'
print '	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
print '	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">'
print '<head><title>HVDRP Serotyper</title>'
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head><body>'

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file!', '')

startS = form["starts"].value
if len(startS) == 0:
	badFile('Error: Specify motif!', '')

truncate = form["truncate"].value

print '<h1> <style="vertical-align:middle">"Purdy" HBV Serotyper Results</h1><br>'

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'Serotyper-' + randomToken + '/'

os.makedirs(climbMini.TEMPFOLDER)

tempF = climbMini.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

S = climbMini.Sequence()
S.load(tempF)
if int(truncate) != 0:
	S.truncate = int(truncate)

Sero = S.seqSerotype(startMotif = startS)

print '<b>Input file:</b> %s<br><br>' % filename
print '<b>S Region Start Motif:</b> %s<br><br>' % startS

print '<div class="dotted_blue" style="border-width:1.5px; color:darkblue; font-size:12px; font-family:monospace; width:66%; margin:2px">'

print 'The "Motif1" column shows the five amino acids (122, 160, 127, 159 and 140) required for serotyping. All five amino acids are not required to determine all serotypes, but all five are shown for each sequence. The "Motif3" column shows the three-letter abbreviations for the amino acid motif. The "Motif Sequence" shows the nucleotide sequence for each of the five amino acids for reference (positions 366-368, 480-482, 381-383, 477-479 and 420-422). Serotypes which cannot be determined are shown as "Unknown". A serotype of "(ad)" or "(ay)" indicates that only the <i>d</i> or <i>y</i> determinant class could be deduced. A summary table of the amino acids required for each serotype appears below the results.<br><br>The serotype is determined as per Figure 3 in Purdy, M. A., Talekar, G., Swenson, P., Araujo, A. and Fields, H. 2007. A new algorithm for deduction of hepatitis B virus surface antigen subtype determinants from the amino acid sequence. <i>Intervirology</i> <b>50</b>: 45-51. DOI: 10.1159/000096312.<br></div><br>'

print '<table border="1" padding="1">'
print '<tr><td><b>Sequence ID</b></td><td><b>Serotype</b></td><td><b>Motif1</b></td><td><b>Motif3</b></td><td><b>Motif Sequence</b></td></tr>'

for i in Sero:
	print '<tr style="background-color:%s;">' % climbMini.HIGHLIGHTSEROTYPES[i[1]]
	for j in i:
		print '<td>%s</td>' % j
	print '</tr>'
print '</table><br>'

print '<table border="1" style="font-family:monospace">'
print '<tr><td align="center"></td><td align="center">122</td><td align="center">160</td><td align="center">127</td><td align="center">159</td><td align="center">140</td><tr>'
print '<tr><td align="center">adr </td><td align="center">K</td><td align="center">R</td><td align="center"></td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">adw2</td><td align="center">K</td><td align="center">K</td><td align="center">P</td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">adw3</td><td align="center">K</td><td align="center">K</td><td align="center">T</td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">adw4</td><td align="center">K</td><td align="center">K</td><td align="center">I or L</td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">ayr </td><td align="center">R</td><td align="center">R</td><td align="center"></td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">ayw3</td><td align="center">R</td><td align="center">K</td><td align="center">T</td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">ayw4</td><td align="center">R</td><td align="center">K</td><td align="center">I or L</td><td align="center"></td><td align="center"></td><tr>'
print '<tr><td align="center">ayw1</td><td align="center">R</td><td align="center">K</td><td align="center">P</td><td align="center">A</td><td align="center"></td><tr>'
print '<tr><td align="center">ayw2</td><td align="center">R</td><td align="center">K</td><td align="center">P</td><td align="center">Not A</td><td align="center">Not S</td><tr>'
print '<tr><td align="center">ayw4</td><td align="center">R</td><td align="center">K</td><td align="center">P</td><td align="center">Not A</td><td align="center">S</td><tr>'
print '</table>'

S.unload()

closePage()

