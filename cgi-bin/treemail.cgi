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

# 15 May 2011

import datetime

starttime = datetime.datetime.now()

import climbMini, cgi, os, re, sys, socket, shutil, subprocess

form = cgi.FieldStorage()

version = "<h1>Pipeline: TreeMail</h1>"

def closePage():
	print '<hr>'
	print '<tt><p align=right>Version 0.2 Beta (May 2011)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'Pipeline: TreeMail\t' )
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
print '<head><title>Pipeline: TreeMail</title>'
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head><body>'

cgi.sys.stdout.flush()

email = form["emailaddress"].value
if len(email) <= 6 or email.find('@') == -1 or email.find('.') == -1:
	badFile("Bad email address!", "")

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file(s)!', '')

filter = form["filter"].value

if filter == "01":
	mode = "TreeMail"
else:
	mode = "TreeMailBootstrap"

print 'Setting up<br>'
cgi.sys.stdout.flush()

randomToken = climbMini.randomStamp()
climbMini.TEMPFOLDER = '/var/www/tmp/'
climbMini.TEMPFOLDER += mode + '-' + randomToken + '/'

os.makedirs(climbMini.TEMPFOLDER)

shutil.copy(mode + ".sh", climbMini.TEMPFOLDER + mode + ".sh")

tempF = climbMini.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

print 'Running script ... please wait ...<br>'
cgi.sys.stdout.flush()

subprocess.call([climbMini.TEMPFOLDER+mode+".sh", filename], stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, cwd=climbMini.TEMPFOLDER)

import os
if not os.path.isfile(climbMini.TEMPFOLDER + 'result.tre'):
	print "Error in input file. Check format."
	closePage()
	sys.exit()

print 'Sending email<br>'
cgi.sys.stdout.flush()

import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders

gmail_user = "user@domain.com"
gmail_pwd = "password"

def mail(to, subject, text, attach1, attach2):
   msg = MIMEMultipart()

   msg['From'] = gmail_user
   msg['To'] = to
   msg['Subject'] = subject

   msg.attach(MIMEText(text))

   part = MIMEBase('application', 'octet-stream')
   part.set_payload(open(attach1, 'rb').read())
   Encoders.encode_base64(part)
   part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(attach1))
   msg.attach(part)

   part = MIMEBase('application', 'octet-stream')
   part.set_payload(open(attach2, 'rb').read())
   Encoders.encode_base64(part)
   part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(attach2))
   msg.attach(part)


   mailServer = smtplib.SMTP("smtp.gmail.com", 587)
   mailServer.ehlo()
   mailServer.starttls()
   mailServer.ehlo()
   mailServer.login(gmail_user, gmail_pwd)
   mailServer.sendmail(gmail_user, to, msg.as_string())
   # Should be mailServer.quit(), but that crashes...
   mailServer.close()

if filter == "02":
	note = "(consensus from bootstrap) "
else:
	note = ""

mail(email, "TreeMail Result: "+filename, "Tree file "+note+ "from input file "+filename+" attached.", climbMini.TEMPFOLDER + "result.tre", climbMini.TEMPFOLDER + "result.txt")

print 'Done<br>'

closePage()

