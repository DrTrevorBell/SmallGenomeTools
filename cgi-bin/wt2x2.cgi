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

version = "Wild-type 2x2"

def closePage():
	print '<br><a href="/wt2x2/index.html">Submit another</a>'
	print '<hr>'
	print '<tt><p align=right>Version 1.0 (February 2012)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'wt2x2\t' )
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
print 'td { text-align:center; }'
print 'td.indivheading { background:lightgreen; text-align:center; color:black }'
print 'td.heading1 { background:lightblue; text-align:center; }'
print 'td.heading2 { background:#ffff99; text-align:center; }'
print 'td.info { text-align:left; }'
print '</style>'
print '</head>'
print '<body>'

infile = form["infile"].file
filename = form["infile"].filename

if (len(filename) == 0):
	badFile('Error: Specify input file.', '')

group1Name = form['name1'].value
group2Name = form['name2'].value

theThreshold = float(form['thresh'].value)
theOffset = int(form['offset'].value)

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/wt2x2-' + randomToken + '/'
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

# process grouping
group1 = form['g1'].value
group2 = form['g2'].value

def makeOtherList(onlyList):
	otherList = ''		# create this list as a string so that it can be processed in the same way as if it were entered

	gg = [climbMini.parseList(onlyList)]
	g = []
        for i in gg:
		tt = []
                for j in i:
			for k in range(j[0], j[1]+1):
				tt.append(k)

	for i in range(S.seqCount()):
		if i+1 not in tt:
			otherList = otherList + str(i+1) + ','

	return otherList[:-1]	# remove trailing ,
# ----------

if len(group2) == 0:
	group2 = makeOtherList(group1)

if len(group1) == 0:
	group1 = makeOtherList(group2)

theGrouping = [group1, group2]

results = S.wt2x2(theGrouping, theThreshold)

if theOffset != 0:
	for i in results:
		i[0] += theOffset - 1

print '<h1>%s</h1>' % version
print '<h2>Results</h2>'

print '<table>'
print '<tr><td class="info">Filename</td><td class="info">%s</td></tr>' % (filename)
print '<tr><td class="info">Group 1 (%s)</td><td class="info">%s</td></tr>' % (group1Name, theGrouping[0])
print '<tr><td class="info">Group 2 (%s)</td><td class="info">%s</td></tr>' % (group2Name, theGrouping[1])
print '<tr><td class="info">Threshold</td><td class="info">%4.3f</td></tr>' % (theThreshold)
print '</table><br>'

print '<table border="1"><tr><td class="heading1">Position</td><td class="heading1">Wild-type</td><td class="heading1">pFET</td><td class="heading1">Chi</td><td class="heading1">Mutants Group 1 (%s)</td><td class="heading1">Mutants Group 2 (%s)</td></tr>' % (group1Name, group2Name)
for i in results:
	print '<tr>'
	print '<td>%i</td>' % i[0]
	print '<td>%s</td>' % i[1]
	print '<td>%4.3f</td>' % i[2]
	if i[3] == 99:
		chi = '-'
	else:
		chi = '%4.3f' % (float(i[3]))
	print '<td>%s</td>' % chi
	print '<td>%s</td>' % i[6]
	print '<td>%s</td>' % i[7]
	print '</tr>'
print '</table>'
for i in results:
	print '<br><table border="1">'
	print '<tr><td colspan="4" class="indivheading">%i</td></tr>' % (i[0])
	print '<tr><td class="heading2">Group</td><td class="heading2">Wild-type</td><td class="heading2">Mutant</td><td class="heading2">Total</td></tr>'
	print '<tr><td>1 (%s)</td><td>%i</td><td>%i</td><td>%i</td></tr>' % (group1Name, i[4][0], i[4][1], i[4][0] + i[4][1])
	print '<tr><td>2 (%s)</td><td>%i</td><td>%i</td><td>%i</td></tr>' % (group2Name, i[5][0], i[5][1], i[5][0] + i[5][1])
	print '<tr><td class="heading2">Total</td><td class="heading2">%i</td><td class="heading2">%i</td><td class="heading2">%i</td></tr>' % (i[4][0] + i[5][0], i[4][1] + i[5][1], i[4][0] + i[4][1] + i[5][0] + i[5][1])

	print '</table>'

closePage()

