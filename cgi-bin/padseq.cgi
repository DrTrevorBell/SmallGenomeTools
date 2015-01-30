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

version = "HBV PadSeq"

def closePage():
	print '<br><a href="/padseq/index.html">Submit another</a>'
	print '<hr>'
	print '<tt><p align=right>Version 1.0 (February 2012)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'PadSeq\t' )
	f.write( 'Request from %s\t' % (cgi.os.environ['REMOTE_ADDR']) )
	f.write( 'Served by %s at %s\t' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0]) )
	f.write( 'Run completed %s\t' % (endtime.strftime("%Y-%m-%d %H:%M")) )
	f.write( '%03.6f seconds\t' % ((endtime - starttime).microseconds / float(1000000)) )
	f.write( '\n')
	f.close()

def badFile(errorText, explainText):
	print "<h1>%s</h1>" % (version)
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

infile1 = form["infile1"].file
filename1 = form["infile1"].filename

infile2 = form["infile2"].file
filename2 = form["infile2"].filename

if (len(filename1) == 0) or (len(filename2) == 0):
	badFile('Error: Specify all input files.', '')

Length = int(form['length'].value)
Genotype = 'Genotype' + form['genotype'].value
Backbone = form['backbone'].value
fragPositions = []
fragPositions.append(int(form['file1pos'].value))
fragPositions.append(int(form['file2pos'].value))

if (fragPositions[0] > Length) or (fragPositions[1] > Length) or (fragPositions[0] <= 0) or (fragPositions[1] <= 0):
	badFile('Error: Bad position/s. Check and resubmit.', '')

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/PadSeq-' + randomToken + '/'
os.makedirs(climbMini.TEMPFOLDER)
tempF = climbMini.TEMPFOLDER + filename1
temp = open(tempF, 'wb')
for line in infile1:
	temp.write(line)
temp.close()
infile1.close()

S = climbMini.Sequence()
S.load(tempF)
S.seqCase()

tempF = climbMini.TEMPFOLDER + filename2
temp = open(tempF, 'wb')
for line in infile2:
	temp.write(line)
temp.close()
infile2.close()

T = climbMini.Sequence()
T.load(tempF)
T.seqCase()

if S.seqCount() != T.seqCount():
	badFile('Error: Number of sequences does not match.', '')

outFilename = '/var/www/tmp/PadSeq-' + randomToken + '/PadSeq.fasta'
outFilenameLink = '/tmp/PadSeq-' + randomToken + '/PadSeq.fasta'
outFile = open(outFilename, 'wb')

def placeFrag(ss, pp, bb):	# accepts the ACTUAL position (1-indexed) as this is what is specified by the user and passed in for non-wraps positions
	return bb[:pp-1] + ss + bb[pp+len(ss)-1:]

for i in range(S.seqCount()):
	backbone = Backbone * Length
	frags = []	# get both fragments into one list for processing later
	frags.append(S.seq[i]['seq'])
	frags.append(T.seq[i]['seq'])
	for j in range(len(frags)):
		if fragPositions[j] + len(frags[j]) -1 > Length:
			backbone = placeFrag(frags[j][:Length-fragPositions[j]+1], fragPositions[j], backbone)
			backbone = placeFrag(frags[j][Length-fragPositions[j]+1:], 1, backbone)		# passing in the ACTUAL position (1-indexed)
		else:
			backbone = placeFrag(frags[j], fragPositions[j], backbone)

	outFile.write(">%s\n" % (S.seq[i]['id'] + "__" + T.seq[i]['id']))
	outFile.write("%s\n" % (backbone))

outFile.close()

print 'Padding completed.<br><br><a href="%s">Download</a> padded FASTA file of all sequences (right-click and "Save As").<br>' % (outFilenameLink)

closePage()

