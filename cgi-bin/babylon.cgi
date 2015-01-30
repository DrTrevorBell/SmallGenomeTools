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

version = "HBV Babylon Translator"

def closePage():
	print '<br><a href="/babylon/index.html">Submit another</a>'
	print '<hr>'
	print '<tt><p align=right>Version 1.0 (January 2012)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'Babylon Translator\t' )
	f.write( 'Request from %s\t' % (cgi.os.environ['REMOTE_ADDR']) )
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
print '<head><title>HVDR %s</title>' % (version)
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '<style type="text/css">'
print 'td.ok { background:green; text-align:center; }'
print 'td.bad { background:red; text-align:center; }'
print '</style>'
print '</head>'
print '<body>'

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file!', '')

proteinName = ['Pre-Core', 'Core', 'Pol' ,'PreS1', 'PreS2' ,'Surface', 'X']
proteinList = []
proteinCount = 0

Length = form['length'].value
Translate = 'translate' in form
Genotype = 'Genotype' + form['genotype'].value
Subs = 'subs' in form

proteinPositions = []
for i in range(1,8):
	if 'include' + str(i) in form:	# has_key
		proteinList.append(proteinName[i-1])
		proteinCount += 1
		proteinPositions.append(i)


coords = []
for i in proteinPositions:
	coords.append((form['s'+str(i)].value, form['e'+str(i)].value))

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'Babylon' + randomToken + '/'
os.makedirs(climbMini.TEMPFOLDER)
tempF = climbMini.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

S = climbMini.Sequence()
S.load(tempF)

S.seqCase()

if Subs:
	for j in S.seq:
		j['seq'] = j['seq'].replace('-', 'N')
		j['seq'] = j['seq'].replace('?', 'N')

print '<b>Input file:</b> %s<br><br>' % filename

checkStart = []
checkEnd = []
for i in range(proteinCount):
	checkStart.append(S.extract(coords[i][0]+'-'+str(int(coords[i][0])+2)))
	checkEnd.append(S.extract(str(int(coords[i][1])-2)+'-'+coords[i][1]))

print '<table border="1">'
print '<tr><td>Sample</td>'
for i in proteinList:
	print '<td style="text-align:center;" colspan="2">%s</td>' % i
print '<td rowspan="2" style="text-align:center;">Length</td></tr><td>Protein</td>'
for i in proteinList:
	print '<td style="text-align:center;">Start</td><td  style="text-align:center;">End</td>'
print '</tr>'

for i in range(S.seqCount()):
	print '<tr><td>%s</td>' % checkStart[0][i][0]
	for j in range(proteinCount):
		out = checkStart[j][i][1]
		if out == 'ATG':
			outClass = 'ok'
		else:
			outClass= 'bad'
		print '<td class="%s">%s</td>' % (outClass, out)

		out = climbMini.robustTranslate(checkEnd[j][i][1])[0]
		if out == '*':
			outClass = 'ok'
		else:
			outClass = 'bad'
		print '<td class="%s">%s</td>' % (outClass, out)

		SL = S.seqLength()[i][1]
		if str(SL) == Length:
			outClass = 'ok'
		else:
			outClass = 'bad'
	print '<td class="%s">%s</td>' % (outClass, str(SL))

	print '</tr>'
print '</table>'

fileList = []
for i in range(proteinCount):
	start = int(coords[i][0])
	end = int(coords[i][1])
	if end < start:
		out = coords[i][0]+'-'+Length+',1-'+coords[i][1]
	else:
		out = coords[i][0]+'-'+coords[i][1]
	outDetail = S.extract(out, aaout=Translate)
	outFilename = Genotype+'-'+proteinList[i]+'.fasta'
	fileList.append(outFilename)
	outFile = open(climbMini.TEMPFOLDER+outFilename, 'w')
	for j in range(S.seqCount()):
		outFile.write('>%s\n%s\n' % (outDetail[j][0]+'-'+proteinList[i], outDetail[j][1]))
	outFile.close()

import zipfile, os
os.chdir(climbMini.TEMPFOLDER)
Z = zipfile.ZipFile(climbMini.TEMPFOLDER+'Babylon.zip', 'w')
for zipName in fileList:
	Z.write(zipName)
Z.close()
S.unload()

for i in range(len(fileList)):
	filename = "/tmp/Babylon" + randomToken + "/" + fileList[i]
	print '<br><a href="%s">Download</a> %s FASTA file<br>' % (filename, proteinList[i])

filename = "/tmp/Babylon" + randomToken + "/Babylon.zip"
print '<br><a href="%s">Download</a> all of the above in one ZIP file<br>' % (filename)

print '<br>The download files  will expire after one hour.<br>'

closePage()

