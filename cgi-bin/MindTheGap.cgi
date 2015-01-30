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

import climbMini
import sys
import cgi
import os
import socket

form = cgi.FieldStorage()

version = "<h1>Mind the Gap</h1>"

def closePage():
	print headerLinks
	print '<hr><input type="button" value="Go Back" onclick="goBack()" /><br>'
	print '<hr>'
	print '<tt><p align=right>Version 0.1 (November 2014)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'MindTheGap\t' )
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
print '<head><title>HVDR Mind the Gap</title>'
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head>'
print '<script type="text/javascript">'

print 'function goBack()'
print '  {'
print '  window.history.back()'
print '  }'
print '</script>'

print '<style type="text/css" rel="stylesheet">'
print '<!--'
print '.hidden { visibility: hidden; }'
print '.unhidden { visibility: visible; }'
print 'tr.main {font-size: 12px; text-align: center; font-family: monospace; }'
print 'td.Correct {background: brown; font-weight: bold; }'
print 'td.Check {background: coral; }'

print 'td.header {text-align: center; background: lightblue; font-weight:bold; font-family: monospace; }'
print 'td.summary1 {text-align: center; background: firebrick; font-family: monospace; }'
print 'td.summary2 {text-align: center; background: chocolate; font-family: monospace; }'
print 'td.ignore {text-align: center; background: black; color:white; font-family: monospace; }'

print 'td.Center {text-align: center; }'
print 'td.o {font-size: 12px; color: black; text-align:center; font-family: monospace; }'

print 'table.fixed {table-layout:fixed; }'
print '-->'
print '</style>'

print '<body>'

headerLinks = ''

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file!', '')

threshold = form['threshold'].value
threshold = threshold.strip()
threshold = threshold.upper()
threshold = float(threshold)/100

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += 'MindTheGap-' + randomToken + '/'
os.makedirs(climbMini.TEMPFOLDER)
os.chdir(climbMini.TEMPFOLDER)

tempF = climbMini.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

S = climbMini.Sequence()
S.load(tempF)
S.seqCase()

d = S.baseDistribution('1-%s' % str(S.originalLength))

print 'Using threshold %5.2f<br><br>' % threshold

# ===== Output gap distribution for all columns with at least one gap =====
for i in range(len(d)):
        if d[i][climbMini.GAP] != 0:
                print 'Position %04i: Gaps: %04i/%04i (%05.2f%%)   Non-gaps: %04i/%04i (%05.2f%%)<br>' % (i+1, d[i][climbMini.GAP], S.seqCount(), float(d[i][climbMini.GAP])/S.seqCount()*100, S.seqCount()-d[i][climbMini.GAP], S.seqCount(), float(S.seqCount()-d[i][climbMini.GAP])/S.seqCount()*100)

print

# ===== Remove sequences with a non-gap in columns with >= threshold gaps =====
removeListDict = {}
removeList = []
T = climbMini.Sequence()
T.fileLoaded = True
T.fileName = 'NoInsertion-%3.2f.fasta' % (threshold)
noInsertionName = T.fileName
for i in range(len(d)):
        if d[i][climbMini.GAP] / float(S.seqCount()) >= threshold:
                for j in range(S.seqCount()):
                        if S.seq[j]['seq'][i] != climbMini.GAP:
                                if S.seq[j]['id'] not in removeListDict:
                                        T.seqAdd([S.seq[j]])
                                        removeListDict[S.seq[j]['id']] = None

for i in removeListDict:
        removeList.append(i)
print '<br>Removing the following %i sequence(s):<br>' % (len(removeList))
for i in removeList:
	print i, '<br>'

print '<br>Sequence count before removing: %i<br>' % S.seqCount()
T.writeFASTA(filename=T.fileName)
T.unload(override=True)

# remove from the list now, rather than in the loop above
S.seqRemoveByID(removeList)

print 'Sequence count after removing:  %i<br>' % S.seqCount()

d = S.baseDistribution('1-%s' % str(S.originalLength))

if S.seqCount() > 0:	# there may be no sequences left
	S.eliminateGapColumns(gapThreshold=threshold)

insertionName = 'Insertion-%3.2f.fasta' % (threshold)
S.writeFASTA(filename=insertionName)

S.unload(override=True)

print '<br>'
print 'Download FASTA file containing sequences without insertion <a href="%s">here</a> (right-click).<br><br>' % ("/tmp/MindTheGap-" + randomToken + '/' + noInsertionName)
print 'Download FASTA file containing sequences with insertion <a href="%s">here</a> (right-click).<br>' % ("/tmp/MindTheGap-" + randomToken + '/' + insertionName)


closePage()

sys.exit()

