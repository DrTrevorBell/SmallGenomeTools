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

import os
import re
import sys
from Bio import SeqIO
import Bio.Seq

PROGRAM = 'CLIMBEDmini'
VERSION = '0.1.20150115a'

BASEREFERENCE = { 'A'     : 'A', 'C'   : 'C',  'G'  : 'G', 'T'   : 'T', '-'  : '-',
                  'AC'   : 'M', 'AG'  : 'R', 'AT'  : 'W', 'CG'  : 'S', 'CT'  : 'Y', 'GT'  : 'K',
                  'ACG'  : 'V', 'ACT' : 'H', 'AGT' : 'D', 'CGT' : 'B',
                  'ACGT' : 'N'}

BASES = ['A', 'C', 'G', 'T']

NONBASES = ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N', '-']

BASECOMPLEMENT = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', '-' : '-', 'N' : 'N',
                  'K' : 'M', 'M' : 'K', 'R' : 'Y', 'Y' : 'R', 'W' : 'W', 'S' : 'S',
                  'B' : 'V', 'V' : 'B', 'D' : 'H', 'H' : 'D', 'X' : 'X',
                  'a' : 't', 't' : 'a', 'c' : 'g', 'g' : 'c',            'n' : 'n',
                  'k' : 'm', 'm' : 'k', 'r' : 'y', 'y' : 'r', 'w' : 'w', 's' : 's',
	          'b' : 'v', 'v' : 'b', 'd' : 'h', 'h' : 'd', 'x' : 'x'}

GAP = '-'

DISAMBIGUATE = { 'M' : 'AC',  'R' : 'AG',  'W' : 'AT',  'S' : 'CG',  'Y' : 'CT',  'K' : 'GT',
		 'V' : 'ACG', 'H' : 'ACT', 'D' : 'AGT', 'B' : 'CGT', 'N' : 'ACGT' }

QUALITYTHRESHOLD = 20
BASESTHRESHOLD = 5
GAPCOLUMNTHRESHOLD = 0.80

AMINOACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

HIGHLIGHTAMINOACIDS = {
'Met': 'green',
'Ter': 'red',
'Xaa': 'gray',
'Xaa': 'gray',
'Ala': '#BFBFFE',
'Arg': '#FF80FE',
'Asn': '#8080FE',
'Asp': '#BFFFFE',
'Cys': '#80FFFE',
'Gln': '#9999FF',
'Glu': '#7A7ACC',
'Gly': '#FEFF80',
'His': '#B2B300',
'Ile': '#FEBFBF',
'Leu': '#CCCC00',
'Lys': '#BFFEBF',
'Phe': '#80FE80',
'Pro': '#FE8080',
'Ser': '#FF66CC',
'Thr': '#00B3B2',
'Trp': '#8FFEBF',
'Tyr': '#FEFFBF',
'Val': '#80FE80'
}

SEROPATTERNS = { 'KR...' : 'adr', 'KKP..' : 'adw2', 'KKT..' : 'adw3', 'KK[I|L]..' : 'adw4', 'RR...' : 'ayr',
		                 'RKT..' : 'ayw3', 'RK[I|L]..' : 'ayw4', 'RKPA.' : 'ayw1', 'RKP[^A][^S]' : 'ayw2', 'RKP[^A]S' : 'ayw4' }

SEROGROUP = { 'K' : '(ad)', 'R' : '(ay)' }

GENOTYPES = {'adw2CGTCA[T|Y|C]...C[A|C]' : 'A1',
             'ayw1CGTCA[T|Y|C]...C[A|C]' : 'A1',
             'adw2CGGCAC...CG' : 'A2',
             'ayw1CGTCA[T|Y|C]...CG' : 'A3',
             'ayw1CG.......CG' : 'A4',
             'adw.TT.......T.' : 'B1/B2/B3/B6',
             'ayw.TT.......T.' : 'B3/B4/B5',
             'adrTT.......C.' : 'C1/C2',
             'ayrTT.......C.' : 'C1/C2',
             'adrCG.......C.' : 'C3',
             'ayrCG.......C.' : 'C3',
             'ayw.TT.......T.' : 'C4',
             'adw2CG.......T.' : 'C5',
             'ayw[^4]CG.......T.' : 'D',
             'ayw4CG.......T.' : 'E',
             'adw4TT.......T.' : 'F1/F4',
             'adw4TT.......C.' : 'F2/F3/H',
             'adw2CG....TAAT.' : 'G' }

HIGHLIGHTSEROTYPES = {
'(ad)'    : '#FEFF80',
'(ay)'    : '#7A7ACC',
'Unknown' : 'gray',
'adr'     : '#FEFFBF',
'adw2'    : '#80FE80',
'adw3'    : '#B2B300',
'adw4'    : '#BFBFFE',
'ayr'     : '#FF80FE',
'ayw1'    : '#80FFFE',
'ayw2'    : '#FF66CC',
'ayw3'    : '#9999FF',
'ayw4'    : '#BFFFFE'
}

S_START = '[A|a][T|t|C|c|G|g][G|g][G|g][A|a][G|A|C|S|g|a|c|s][A|G|R|a|g|r][A|G|R|a|g|r][C|c][A|a][C|T|c|t]'
B_START = '[A|a|C|c|T|t][T|t|C|c][G|g|A|a][C|T|c|t][A|a][A|a][C|c][T|t][T|t][T|t][T|t][T|t][C|c][A|a][C|c][C|c][T|t][C|c][T|t][G|g]'
TEMPFOLDER = '/tmp/'
MAX_FRAGMENTS = 12

if __name__ == '__main__':
	interactive = True
else:
	interactive = False

def motd():
	divider = '-' * 40 + '\n'
	return divider + 'This is %s %s\n\nDependencies:\n\tBioPython from http://www.biopython.org\n\tABIFReader.py from http://www.interactive-biosoftware.com/open-source/ABIFReader.py\n' % (PROGRAM, VERSION) + divider
# ---------------

def error(message):
	out = 'Error [%s]: %s' % (sys._getframe(1).f_code.co_name, message)

	if interactive:
		sys.stderr.write(out + '\n')
	else:
		sys.exit(out)	# generates errorlevel of 1 (check with 'echo $?' after running); no tracebacks	# was 'raise Exception(out)'

def warning(message):
	# Name of calling method: http://bytes.com/topic/python/answers/665113-how-can-i-know-name-caller
	out = 'Warning [%s]: %s ' % (sys._getframe(1).f_code.co_name, message)
	if interactive:
		print (out)
	else:
		sys.stderr.write(out)
# ---------------

def parseList(tempcc, aa=False, mappingPos=1):	# default mapping is 1:1
	'''Returns a start and end value for each distinct entity'''
	# For example: 1800-1810 returns 1800 and 1810
	# 1820 returns 1820 and 1820
	ll = tempcc.split(',')
	tempNucList = []
	for i in ll:
		rr = i.find('-')
		if rr != -1:
			tt = i.split('-')
			start = tt[0]
			end = tt[1]
		else:
			start = end = i

		try:
			startPos = int(start)
			endPos = int(end)
			if aa:
				startPos = startPos * 3 - 2		# start position
				endPos = endPos * 3			# do not subtract 2 -- triplets
			tempNucList.append([startPos - mappingPos + 1, endPos - mappingPos + 1])	# a list of the positions requested
		except:
			error('%s: Bad nucleotide positions' % tempcc)
			return False

	return tempNucList
# ---------------

def randomStamp():
	import datetime, random
	return "%s%02i" % (datetime.datetime.strftime(datetime.datetime.now(),'%Y%m%d%H%M%S'), random.random()*999)
# ----------------

def contig(F, R):
	'''Generate a contig from two input sequences'''
	# Method requires two sequences objects which may optionally have been trimmed; objects required because quality scores are required
	# Processes the first sequence in each Sequence object
	# Ancestry: stand-alone qualTrim.py modified into stand-alone Contig.py

	# Tested with SHH061A BCP forward and reverse sequences; after _5_10 trimming and needle alignment, base T at 220 with QS of 54 was aligned with gap

	# This needs a forward and reverse file, called F and R -- must parameterize this
	# Pass in the entire object; only the first sequence of the object used

	from Bio.Emboss.Applications import NeedleCommandline
	from Bio import AlignIO
	import subprocess, sys

	R.seqRevComp()

	randomStampToken = randomStamp()	# same token for forward and reverse files

	tempOutF = TEMPFOLDER + 'tempOutF-' + randomStampToken
	tempOutR = TEMPFOLDER + 'tempOutR-' + randomStampToken

	F.writeFASTA(tempOutF)
	R.writeFASTA(tempOutR)

	# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc74
	# Using BioPython tools to call needle with prepared files
	# Uses subprocess to pipe stdout, which avoids making a temporary output file

	cline = NeedleCommandline(asequence=tempOutF, bsequence=tempOutR, gapopen=10, gapextend=0.5, stdout=True, auto=True)
	child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
	align = AlignIO.read(child.stdout, "emboss")

	# Quality scores always preceed bases
	pos1 = -1	# track the actual base position regardless of gaps
	pos2 = -1
	gap = '-'
	cons = ''

	# use a moving window around the current pair of bases
	# the quality of the bases in this window will be conisidered
	# SHH009A reverse sequence contains an an unpaired Y base with quality of 5, but this base is valid and should not be deleted
	# arbitrarily deleting unpaired bases is not advisable; possibly lower thresholds on non-ACGT bases

	CW = BASESTHRESHOLD 		# CW ~ ContigWindow
	QT = QUALITYTHRESHOLD		# shorter to type and only one change here if module name is changed

	pairBaseGood = 1.0
	pairBaseBase = 0.75
	pairBaseNon = 0.67
	pairGapBase = 0.50
	pairBaseBasePoor = 0.33
	PairNon = 0.25

	CQ = 0.0		# Contig Quality

	for i in range(len(align[0].seq)):
		base1 = align[0].seq[i]
		base2 = align[1].seq[i]

		if base1 != gap:
			pos1 += 1
		if base2 != gap:
			pos2 += 1

		qual1 = F.QS[pos1 + F.trimLeft]
		qual2 = R.QS[pos2 + R.trimLeft]

		window1 = [pos1 + F.trimLeft - CW, pos1 + F.trimLeft + CW]
		window2 = [pos2 + R.trimLeft - CW, pos2 + R.trimLeft + CW]

		if window1[0] < 0:
			window1[0] = 0
		if window1[1] > len(F.seq[0]['seq']):
			window1[1] = len(F.seq[0]['seq'])

		if window2[0] < 0:
			window2[0] = 0
		if window2[1] > len(R.seq[0]['seq']):
			window2[1] = len(R.seq[0]['seq'])

		baseWin1 = F.seq[0]['seq'][window1[0]:window1[1]+1]
		baseWin2 = R.seq[0]['seq'][window2[0]:window2[1]+1]

		qualWin1 = sum(F.QS[window1[0]:window1[1]+1]) // (CW * 2 + 1)	# integer division
		qualWin2 = sum(R.QS[window1[0]:window1[1]+1]) // (CW * 2 + 1)

		# print "pos1:%03i qual1:%i base1:%s window1:%03i-%03i basewin1:%s qualwin1:%i | pos2:%03i qual2:%i base2:%s window2:%03i-%03i basewin2:%s qualwin2:%i" % (pos1, qual1, base1, window1[0], window1[1], baseWin1, qualWin1, pos2, qual2, base2,  window2[0], window2[1], baseWin2, qualWin2)

		# using (quality score of base (qual) OR quality score of window (window) means that
		# an indiviual good quality base is never lost, but also that one poor quality base
		# surrounded by good quality bases is not removed

		# base1 is a gap and base2 to is not a gap and (base2 quality score OR quality score of window2 >= quality threshold)
		if base1 == gap and base2 != gap and (qual2 >= QT or qualWin2 >= QT):
			# print i, base1, qual2, base2
			cons += base2
			CQ += pairGapBase
		# base2 is a gap and base1 to is not a gap and (base1 quality score OR quality score of window1 >= quality threshold)
		elif base1 != gap and base2 == gap and (qual1 >= QT or qualWin1 >= QT):
			# print i, qual1, base1, base2
			cons += base1
			CQ += pairGapBase
			# if one of the bases is a gap and the one which is not a gap has BOTH a quality score and a quality score of window below threshold, then nothing is added to the consensus
		# Bases do not match:
		elif base1 != base2 and base1 != gap and base2 != gap:		# extra clauses just to ensure that neither are gaps
			# If both are ACGT then use the one with the highest quality score as long as it is above qualityThreshold
			# If only one is ACGT, then use that one (currently regardless of quality)
			# If neither are ACGT, then use the one with the best quality (above quality threshold)
			if base1 in BASES and base2 in BASES:
				if qual1 >= qual2 and qual1 >= QT:
					cons += base1
					CQ += pairBaseBase
				elif qual2 >= QT:
					cons += base2
					CQ += pairBaseBase
				# both are ACGT but neither are above qualityThreshold
				elif qual1 > qual2:	# base1 quality is better
					cons += base1
					CQ += pairBaseBasePoor
				elif qual2 < qual1:	# base2 quality is better
					cons += base2
					CQ += pairBaseBasePoor
				elif qualWin1 > qualWin2:	# base1 window quality is better; will only execute if bases are = quality
					cons += base1
					CQ += pairBaseBasePoor
				else:
					cons += base2
					CQ += pairBaseBasePoor
			elif base1 in BASES and qual1 >= QT:	# base1 is ACGT so base2 is not
				cons += base1			# use base1, regardless of quality
				CQ += pairBaseNon
			elif base2 in BASES and qual2 >= QT:	# base2 is ACGT so base1 is not
				cons += base2			# use base2, regardless of quality
				CQ += pairBaseNon
			else:
				if qual1 >= QT:			# it's not ACGT, but if quality score is above threshold, add it
					cons += base1
					CQ += pairNon
				elif qual2 >= QT:		# ditto; if this is not met, then nothing is added
					cons += base2
					CQ += pairNon
		else:
			cons += base1		# bases are the same so use either -- erroneously takes (took?) "-" from "-|Y"
			CQ += pairBaseGood
		# print cons

	# print "Consensus Length %i" % (len(cons))

	# c = 0
	# for i in cons:
	#	if i not in BASES:
	#		c += 1
	# print "%i of %i (%3.2f%%) bases are ambiguous" % (c, len(cons), float(c)/len(cons) * 100)

	# return contig, CW, aligned forward and reverse sequences, original forward and reverse sequences
	# the original, untrimmed sequence is not available, as this was trimmed when the chromatogram was loaded

	return [[">Contig_%s_%s" % (F.seq[0]['id'][1:], R.seq[0]['id'][1:]), cons],  CQ/len(cons), [">Aligned_Forward_%s" % (F.seq[0]['id'][1:]), str(align[0].seq)], [">Aligned_Reverse_%s]" % (F.seq[0]['id'][1:]), str(align[1].seq)], [F.seq[0]['id'], F.seq[0]['seq']], [R.seq[0]['id'], R.seq[0]['seq']]]

	# Also highlight bases with low quality scores for reference
	# Cannot simply remove base:- pairing arbitrarily if it is found in a long sequence of gaps (as at the start and end of the sequence)
	# Can do it if there position is surrounded by a certain number of base:base pairings -- a threshold
	# Will work well then in the cases examined where base:- should be the base (base with high quality score) or should be removed (base with low quality score)

	# MUST ALSO ADD QUALTRIM PRETTY OUTPUT
	# Output stats: F, R, Contig: As, Gs, Ys, Gaps, Length, etc.
	# Final output must (be/include) a Sequence object?
# ---------------

def robustTranslate(ss):
	'''Use BioPython functionality to translate codons into amino acids'''
	# This method returns a tuple containing the one- and three-letter amino acid codes for the input codon
	# However, this method does not break if the codon cannot be translated
	# as is the case when gaps are present, for example
	# The function will return None if the input parameter is not exactly three characters in length
	if len (ss) != 3:
		return None

	import Bio.SeqUtils

	if ss.find('-') >= 0:
		return ('-', '---')	# means that at least one position is a gap

	ss = ss.upper()

	ss = ss.replace('U', 'T')	# so that translation to three-letter code works

	for i in range(3):
		if ss[i] not in BASES + NONBASES:
			return ('#', '###')	# means that at least one position is not a base or an ambiguous base

	single = Bio.Seq.translate(ss)

	return (single, Bio.SeqUtils.seq3(single))
# ---------------

class Sequence:

	def __init__(self):
		self.fileLoaded = False
		self.fileName = ''
		self.changed = False;
		self.seq = []	# dictonary to store description and sequence
		self.fileType = ''
		self.QS = []
		self.originalCalls = ''
		self.trimThreshold = None
		self.trimLeft = None
		self.originalLength = None
		self.trimRight = None
		self.trimGoodBases = None
		self.trimGoodQuality = None
		self.linenumbers = 0
		self.truncate = 20	# always applied; set to large value to 'disable'
		self.totalCount = 0
		self.mergeSlide = False
	# ---------------

	# def __str__(self):
	#	return "OK"
	# ---------------

	def seqLength(self, cc=''):
		'''Return sequence ID and length'''
		cc = cc.split()
		out = []
		# no marker for data output -- easy copy and paste elsewhere
		if not self.fileLoaded:
			error('No file loaded')
		else:
			for i in range(len(self.seq)):
				out.append((self.seq[i]['id'][:self.truncate], len(self.seq[i]['seq'])))
			# ACTIVE command; SORT
		return out
	# ---------------

	def status(self):
		'''Return status'''
		if not self.fileLoaded:
			error('No file loaded')
		else:
			return (self.fileName, len(self.seq), self.fileType, self.trimThreshold, self.trimLeft, self.trimRight, self.changed)
	# ---------------

	def load(self, cc, ABIF=False, autoTrimThreshold=0, minGoodBases=BASESTHRESHOLD, minGoodQuality=QUALITYTHRESHOLD, filterPattern='.+', filterInclude=True):
		'''Load sequence data from a file'''
		# minGoodBases and minGoodQuality are ignored when autoTrimThreshold is specified
		# filterPattern to filter sequence IDs according to regular expression
		# filterInclude = True includes only sequences matching the filter pattern
		# filters are only applicable to FASTA files, as these can contain multiple sequences
		# filters will be ignored when an ABIF file is specified (this may change in future, but the benefit
		# to filting out an ABIF file is obscure -- the file would not be loaded if the filter did not match,
		# which may be the point in that case)
		### cc = cc.split()
		if len(cc) < 1:
			error('Load what?')

		cc = [cc]
		### if len(cc) > 1:
			### error('Specify one file to load')

		# Check reporting of following ...
		if self.fileLoaded:
			error('Cannot load file %s: File %s already loaded' % (cc[0], self.fileName))

		try:
			if cc[0][0] in ["'", '"'] and cc[0][len(cc[0])-1] in ["'", '"']:
				cc[0] = cc[0][1:len(cc[0])-1]
			if ABIF:
				import ABIFReader
				try:
					f = ABIFReader.ABIFReader(cc[0])
				except:
					print "Bad chromatogram file: %s" % os.path.split(cc[0])[1]	# clumsy but effective
					sys.exit(1)
			else:
				f = open(cc[0], 'r')
		except IOError:
			error('Error opening file %s' % cc[0])
		else:
			self.fileLoaded = True
			self.fileName = cc[0]

			if not ABIF:
				for record in SeqIO.parse(f, "fasta"):
					# id repeated as first part of description, so storing description only
					# length no longer stored here; can be determined on the fly
					# filter each input sequence ID
					self.totalCount += 1	# total all sequences in the file to output include/exclude if relevant
					filterFound = (re.search(filterPattern, record.description) != None)	# return True of found
					if filterFound == filterInclude:
						self.seq.append({'id': record.description, 'seq': str(record.seq)})
						self.originalLength = len(str(record.seq))
						self.trimLeft = 0
						self.trimRight = 0
						# (filterfound == True and filterInclude == True) or (filterFound == False and filterInclude == False)
						# are the two conditions under which the sequence should be included
						# condition above implements this in one condition
				if len(self.seq) < 1:
					if self.totalCount > 1:
						print "All sequences excluded by filter pattern: %s" % filterPattern
					else:
						print "Bad FASTA file: %s" % os.path.split(cc[0])[1]	# clumsy but effective
					sys.exit(1)
				self.fileType = 'FASTA'
			else:
				tempID = f.getData('SMPL')	# Sample ID used as FASTA identifier
				if len(tempID) == 0 or tempID == None:
					tempID = cc[0]		# Filename instead
				base = f.getData('PBAS')	# Base Calls
				qual = f.getData('PCON')	# Quality Scores
				self.date = f.getData('RUND')	# Start run date
				self.QS = []
				for i in qual:
					self.QS.append(ord(i))
				self.originalCalls = base	# preserve the original base calls (as required by fragmentmerger)

				if autoTrimThreshold > 0:
					# This should really be an average of the first X bases?
					# For example: first base QS = 21, then next 10 bases are below 10 ...
					# Or ignore first base?
					# This will break if autoTrimThreshold is very high -- no bases will be good enough
					Qleft = 0
					while self.QS[Qleft] < autoTrimThreshold:
						Qleft += 1
					Qright = len(base) - 1	# indexed from 0
					while self.QS[Qright] < autoTrimThreshold:
						Qright -= 1
					if Qleft <= 0:		# first base is above threshold
						Qleft = 1
					if Qright >= len(base):
						Qright = len(base) - 2	# ?
					base = base[Qleft-1:Qright+1]
	 				self.trimThreshold = autoTrimThreshold
					self.trimLeft = Qleft
					self.trimRight = Qright

					# print self.QS[0:20]
					# print Qleft, Qright
					# print base
				else:
					# No autotrim -- using minGoodBases and minGoodQuality
					c = 0	# counter
					i = 0	# position in sequence
					while i < len(self.QS) and c < minGoodBases:		# still in the sequence and counter < minimum consecutive bases
						if self.QS[i] >= minGoodQuality:
							c += 1
						else:
							c = 0
						i += 1

					if i >= len(self.QS):		# reached the right of the string
						f.close()
						return None
						# error('No good window found searching from the left')	# this could happen if minGoodQuality is too high
					else:
						self.trimLeft = i - minGoodBases	# start of 'good' sequence is X bases 'back' from i

					c = 0	# counter
					i = len(self.QS)-1			# position in sequence, indexed from 0, starting from the right
					while i > 0 and c < minGoodBases:	# still in the sequence and counter < minimum consecutive bases
						if self.QS[i] >= minGoodQuality:
							c += 1
						else:
							c = 0
						i -= 1

					if i == 0:			# reached the left of the string
						f.close()
						return None
						# error('No good window found searching from the right')	# this could happen if minGoodQuality is too high
					else:
						self.trimRight = len(self.QS) - i - minGoodBases - 1

					self.originalLength = len(base)
					base = base[self.trimLeft:len(self.QS)-self.trimRight]
					self.trimGoodBases = minGoodBases
					self.trimGoodQuality = minGoodQuality
					tempID += '_(Trimmed_%i_%i)' % (minGoodBases,minGoodQuality)	# no space so trimmed annotation part of ID

				self.seq.append( {'id': tempID, 'seq': base}) ## MARK A ##
				self.fileType = 'ABIF'

			f.close()

			return self.status()
	# ---------------

	def save(self, cc, overWrite=False):
		'''Save sequence data to a file'''
		# Saves as FASTA
		# Saving with a new name makes that name the 'active' name
		cc = cc.split()
		if len(cc) < 1:
			error('Save what?')

		if len(cc) > 1:
			error('Specify filename to save')

		if not self.fileLoaded:
			error ('Cannot save file %s: No file loaded' % (cc[0]))

		if cc[0][0] in ["'", '"'] and cc[0][len(cc[0])-1] in ["'", '"']:
			cc[0] = cc[0][1:len(cc[0])-1]

		if os.path.exists(cc[0]) and not overWrite:
			error('File exists; specify overWrite=True to overwrite')

		try:
			f = open(cc[0], 'w')
		except IOError:
			error('Error creating file %s' % cc[0])
		else:
			self.changed = False
			self.fileName = cc[0]

			for i in self.seq:
				if i['id'] is not None:			# do not write out empty entries
					f.write('>' + i['id'] + '\n')	# id and description
					f.write(i['seq'] + '\n')	# sequence data

			f.close()

			return self.status()
	# ---------------

	def unload(self, override=False):
		'''Unload file'''
		if self.fileLoaded:
			if self.changed and not override:
				warning('Data changed; specify override=True to unload')
				return None
			t = self.fileName
			self.__init__()
			return None
		else:
			error('Cannot unload: No file loaded')
	# ---------------

	def extract(self, cc, aain=False, aaout=False, mapping=1, inData=''):
		'''Return proteins from a nucleotide sequence'''
		# Base distributions are not returned by this method, as it would be difficult
		# to encapsulate totals for each base for each position requested,
		# and difficult to disentangle these (map them back to positions);
		# the same may be true of the bases themselves, but these are only single characters
		# and produce a 'pattern' or 'motif' which can be processed
		#
		# 04 March 2011: added inData parameter
		# this must be a list (to be interrated over),
		# so to process only one item, pass in inData=[S.seq[0]]
		# several items can be processed with inData=[S.seq[0], S.seq[1]]
		# always returns a list
		#
		cc = cc.split()
		out = []
		if len(cc) < 1:
			error('Specify nucloetide position/s')

		if len(cc) > 1:
			error('Bad position syntax')

		if cc[0].find(' ') != -1:
			error('Spaces not permitted in positions')

		if aain:
			nucList = parseList(cc[0], mappingPos=mapping, aa=True)
		else:
			nucList = parseList(cc[0], mappingPos=mapping)

		if inData == '':
			inData = self.seq

		for i in inData:
			t = i['id'][:self.truncate]
			tempout = ''
			for j in nucList:
				tempout += i['seq'][j[0]-1:j[1]]

			if aaout:
				out.append((t, Bio.Seq.Seq(tempout).translate().data))
			else:
				out.append((t, tempout))

		return out
	# ---------------

	def nucCopy(self, cc):
		'''Crop sequences'''
		self.changed = True
		c = 0
		for i in self.extract(cc):
			self.seq[c]['seq'] = i[1]	# update sequence data only
			c += 1
	# ---------------

	def find(self, cc, aain=False, aaout=False, readingframe=0):
		'''Return all occurrences of a regular expression in a nucleotide sequence'''
		# Returns nucleotide position and, if relevant, amino acid position
		# find over the gap
		# does not return anything if search string not found in sequence (report)
		# returns multiple hits -- some way to limit?
		# output in format which can be added as a new object?
		# removed "context" in output as this did not work well with aain and aaout...
		# supplying aain=False aaout=True with context makes no sense
		if not self.fileLoaded:
			error('No file loaded')
		if len(cc) < 1:
			error('Specify regular expression to find')

		out = []

		outLength = len(cc)
		if aain != aaout:
			if aain and not aaout:
				outLength = outLength * 3
			else:
				outLength = outLength / 3	# not aain and aaout

		for i in self.seq:
			iiid = i['id']
			if aain:
				# Cannot translate gapped sequences; replace - with N for translation to "X"
				ii = i['seq'][readingframe:].replace('-', 'N')
				ii = Bio.Seq.Seq(ii).translate().data
			else:
				ii = i['seq']
			findTemp = re.finditer(cc, ii)
			for foundItem in findTemp:
				startPos = foundItem.start()
				endPos = foundItem.end()
				if startPos < 0:
					startPos = 0
				if endPos > len(ii):
					endPos = len(ii)
				pos = foundItem.start()+1

				if aain:
					nucpos = pos * 3 - 2 + readingframe
					aapos = pos
				else:
					nucpos = pos
					aapos = pos / 3

				if aain != aaout:
					if aain:
						startPos = pos * 3 - 3 + readingframe	# -3 not -2 ... ?
						endPos = startPos + outLength
						lettersOut = i['seq'][startPos:endPos]
					else:
						lettersOut = Bio.Seq.Seq(ii[startPos:endPos]).translate().data
				else:
					lettersOut = ii[startPos:endPos]

				out.append((iiid[:self.truncate], lettersOut, nucpos, aapos))
		return out
	# ---------------

	def seqSlide(self, SStr, SPos):
		'''Place SStr starting at position SPos'''
		self.changed = True
		c = 0
		for i in self.seq:
			if (SPos < 1) or (SPos > len(i['seq'])):
				error('Invalid position')
			start = i['seq'].find(SStr)
			newStart = SPos - 1 - start   # Python indexes strings from zero
			if start < 0:
				warning('%s not found in Sequence %03i (%s)' % (SStr, c, i['id']))
				# return None	# Continue processing rather ...
			else:
				i['seq'] = i['seq'][-newStart:] + i['seq'][:-newStart]
			c += 1
	# ---------------

	def basePercentage(self, cc, byChar=False, percentage=True, ignoreCase=True):
		'''Return count or percentage of base(s) or fragment (no regular expressions)'''
		# cc = str(cc.split()[0])
		out = []
		for i in self.seq:
			if percentage:
				denom = float(len(i['seq']))	# for division later
			else:
				denom = 1
			if ignoreCase:
				cc = cc.upper()
			if not byChar:
				if ignoreCase:
					ccCount = cc.upper()
					count = i['seq'].upper().count(ccCount)
				else:
					ccCount = cc
					count = i['seq'].count(ccCount)
				out.append((i['id'][:self.truncate], ccCount, count / denom))
			else:
				for j in cc:
					if ignoreCase:
						count = i['seq'].upper().count(j)
					else:
						count = i['seq'].count(j)
					out.append((i['id'][:self.truncate], j, count / denom))
		return out
	# ---------------

	def seqRemoveByIndex(self, cc):
		'''Remove sequence (indexed from zero)'''
		# Remove a range or from a list? Numbering?
		self.changed = True
		if cc < 0 or cc > len(self.seq) - 1:
			error('Out of range')
		self.seq.remove(self.seq[cc])
	# ---------------

	# remove sequences by length, by first/last N, odds, evens ...

	def seqRemoveByID(self, cc, retain=False, verbose=False, escape=False):
		'''Remove sequence with regex match against sequence ID'''
		# cc is a list of the regular expressions to be included or excluded
		# verbose for debugging only
		# escape will escape the search pattern; useful for GenBank ID strings which contain |
		if not isinstance(cc, list):
			return None
		self.changed = True
		count = 0
		i = 0
		found = False			# if item is empty, found needs a value later
		while i < len(self.seq):
			if not retain:
				for item in cc:
					if escape:
						tempItem = re.escape(item)
					else:
						tempItem = item
					found = re.search(tempItem, self.seq[i]['id'])
					if found == None:
						found = False
					else:
						found = True
					if found:
						if verbose:
							print item
						self.seqRemoveByIndex(i)
						count += 1
						i -= 1		# one item removed; incremented outside loop
						break		# no need to check the rest of the list
			else:
				for item in cc:
					if escape:
						tempItem = re.escape(item)
					else:
						tempItem = item
					found = re.search(tempItem, self.seq[i]['id'])
					if found == None:
						found = False
					else:
						found = True
					if found:
						break
				if not found:
					if verbose:
						print item
					self.seqRemoveByIndex(i)
					count += 1
					i -= 1			# one item removed; incremented outside loop

			i += 1
		return count			# return number of entries delete (regardless of value of retain)
	# ---------------

	def seqAdd(self, cc):
		'''Add sequences from another object (.seq)'''
		# Must be added from other objects as lists if NewS is empty: NewS.seqAdd([OldS.seq[0]])
		# Or just NewS.seqAdd(OldS.seq) ... ?
		# Can add results of a search: New.seqAdd (OldSeq.seqFind('GAGGAC', 0))
		# This only adds the found sequence not the entire sequence
		# self.fileloaded is not true -- must run seqLength ...
		self.changed = True
		self.seq.extend(cc)
	# ---------------

	def seqCount(self):
		'''Return the number of sequences'''
		if self.fileLoaded:
			return len(self.seq)
		else:
			error('No file loaded')
	# ---------------

	def seqCase(self, lower=False):
		'''Change case of sequence data'''
		self.changed = True
		for i in self.seq:
			if lower:
				i['seq'] = i['seq'].lower()
			else:
				i['seq'] = i['seq'].upper()
	# ---------------

	def seqSR(self, search, replace):
		'''Replace "search" with "replace" in sequence'''
		self.changed = True
		for i in self.seq:
			i['seq'] = re.sub(search, replace, i['seq'])
	# ---------------

	def seqDegap(self):
		'''Remove gaps by calling seqSR'''
		self.seqSR('-', '')
	# ---------------

	def seqRevComp(self, rev=True, comp=True):
		'''Reverse and/or complement sequences'''
		# Might be useful to be able to do one or the other of these on only SOME sequences?
		# Might be useful to simply RETURN the reversed and/or complemented sequence?
		for i in self.seq:
			if rev:
				i['seq'] = i['seq'][::-1]
				self.trimLeft, self.trimRight = self.trimRight, self.trimLeft	# swap
				self.QS = self.QS[::-1]		# Reverse quality scores
			if comp:
				tt = ''		# strings are immutable
				for j in i['seq']:
					tt += BASECOMPLEMENT[j]
				i['seq'] = tt
	# ---------------

	def outFASTA(self):
		out = ''
		for i in self.seq:
			out += '>' + i['id'] + '\n' + i['seq'] + '\n'
		return out
	# ---------------

	def writeFASTA(self, filename):
		try:
			outFile = open(filename, 'w')
		except IOError:
			error('Cannot create file %s' % filename)
		outFile.write(self.outFASTA())
		outFile.close()
	# ---------------

	def countMotif(self, motifPos, motifSeq, match, motifMapping = 1):
		'''Count mutations -- that is, where mutSeq == sequence data'''
		count = 0
		out = []
		pp = self.extract(motifPos, mapping=motifMapping)
		for i in pp:
			if match:
				if i[1] == motifSeq:
					count += 1
					out.append(i)
			else:
				if i[1] != motifSeq:
					count += 1
					out.append(i)
		return (count, out)

	# ----------------

	def baseDistribution(self, positions, set1=BASES, set2=NONBASES, distributionMapping = 1):
		'''Return base distribution values (not percentages) at a given position'''
		# Accepts a parseList-able list of positions and returns a list of dictionaries
		# Processes vertically over all sequences
		out = []

		for i in parseList(positions):
			for j in range(i[0], i[1]+1):
				c = {}
				for k in set1:
					c[k] = 0
				for k in set2:
					c[k] = 0
				for k in self.seq:
					cc = k['seq'][j - distributionMapping]
					if cc in set1 + set2:
						c[cc] += 1
				out.append(c)

		return out	# returned as a list of dictionary; each position is one dictionary
	# ---------------

	def translateLoci(self, positions, translateMapping = 1):
		'''Returns amino acids from all three reading frames for each position'''
		# (Bio.SeqUtils.seq3 converts a Bio.Seq.translate single-letter amino acid to a three-letter amino acid)
		# (Was unable to find a method to to return the full-name of the amino acid)
		# Loop over all positions and extract first reading frame
		# Repeat for other reading frames; therefore three reading frames output per position
		# List of dictionaries is read for output with requiring transposition (as is required for baseDistribtuon)
		# Translation in all three reading frames for all positions is done for completeness and reference
		# This is strictly not necessary for positions which vary by less than three
		out = []
		# Looping over reading frames, then sequences, then positions
		for i in range(3):	# three reading frames per codon
			for j in self.seq:
				c = {}	# prepare the dictionary for this 'row'
				for k in parseList(positions):	# loop over all positions requested for the ith reading frame
					for l in range(k[0], k[1]+1):
						# Returned as a list of dictionaries (as baseDistribution) to allow for easier (modular) output of the table
						# The above turned out not to be possible because of the increased complexity required to store
						# the output from this method; a reading frame, the sequence ID and codons for each position are required
						# The method mentioned above is returing totals for each base, so the data are summarized
						# This method returns output for each reading frame for each sequence for all positions
						# Therefore, a list of list is returned: the inner listcontains the reading frame (0, 1 or 2), the sequence ID
						# and the dictionary of position:codon pairs
						# The key of the dictionary is the actual position number
						# All three reading frames are reaturned, even though some may be empty
						# Returns three reading frames, where position is the first, second and third base:: X--. -X-, --X
						# Reading frames moves 'backwards', but the actual base at the position in question
						# moves from position 1 to 2 to 3, which is more intuitive

						actualPos = l - translateMapping
						start = actualPos - i	# get start position from i
						end = start + 2	# end position is always 2 positions downstream
						if start >= 0:  # upper bound check not required; lower bound check against 0
							# put the codon into dictionary where key is position
							# The translation into an amino acid can be done by the routine
							# which outputs the data, which must consider that not all codons can be translated
							# Codons containing gaps, for exampel, cannot be translated
							# The actual translation into a three-letter amino acid code is done with:
							# import Bio.SeqUtils;  Bio.SeqUtils.seq3(Bio.Seq.translate('AUG'))
							# The outer 'seq3' method can be omitted if one-letter amino acid codes are required
							# Translating functionality is provided by the class method robustTranslate here
							c[l] = j['seq'][start:end+1]

				out.append([i, j['id'][:self.truncate], c])

		return out
	# ---------------

	def eliminateGapColumns(self, gapThreshold=GAPCOLUMNTHRESHOLD, verbose=False):
		'''Eliminate columns containing >= threshold proportion of GAPs'''
		gapList = []
		thresholdList = []

		startingLength = self.seqLength()[0][1]

		for i in range(self.seqLength()[0][1]):		# aligned sequences; all have the same length
			column = ''
			for j in self.seq:
				column += j['seq'][i]
			gapRatio = column.count(GAP) / float(len(self.seq))
			if gapRatio >= gapThreshold:
					gapList.append(i)
					thresholdList.append(gapRatio)

		eliminatedCount = 0
		for i in gapList:
			for j in self.seq:
				j['seq'] = j['seq'][:i - eliminatedCount] + j['seq'][i + 1 - eliminatedCount:]
			eliminatedCount += 1


		if verbose:
			return ['Gap-columns eliminated: %i (%3.2f%%)' % (len(gapList), (len(gapList)/float(startingLength)*100)),
				'Gap threshold: %3.2f' % (gapThreshold),
				'Start length: %i' % (startingLength),
				'End length: %i' % (self.seqLength()[0][1])]
	# ---------------

	def disambiguateColumns(self, verbose=False):
		'''Replace/correct/disambiguate ambiguous bases in mono-base gap-free columns'''
		monobaseCount = 0
		monobaseColumns = 0

		self.seqCase()	# convert to uppercase

		for i in range(self.seqLength()[0][1]):
			bb = []
			column = ''
			for j in self.seq:
				column += j['seq'][i]
			for j in BASES:
				bb.append(column.count(j) / float(len(self.seq)))	# percentage of column which is each basea

			if column.count(GAP) != 0:	# if a gap is present, exclude this column
				continue

			if (bb[1] + bb[2] + bb[3] == 0) and (bb[0] != 1):	# only A present, but not at 100%
				BB = 'A'
			elif (bb[0] + bb[2] + bb[3] == 0) and (bb[1] != 1):	# only C present, but not at 100%
				BB = 'C'
			elif (bb[0] + bb[1] + bb[3] == 0) and (bb[2] != 1):	# only G present, but not at 100%
				BB = 'G'
			elif (bb[0] + bb[1] + bb[2] == 0) and (bb[3] != 1):	# only T present. but not at 100%
				BB = 'T'
			else:
				BB = None					# column contains more than one base-type

			if BB != None:
				monobaseColumns += 1
				for j in column:
					if j != BB:		# if j is the ambiguous base (j is not equal to the base)
						if BB in DISAMBIGUATE[j]:		# the base at this position is one of the disambiguous ones
							monobaseCount += 1
							# print i, column, '<br>'
							# print column[:column.index(j)] + BB + column[column.index(j) + 1:], '<br>'
							# print j, ' becomes ', BB, '<br>'
							self.seq[column.index(j)]['seq'] = self.seq[column.index(j)]['seq'][:i] + BB + self.seq[column.index(j)]['seq'][i + 1:]
							# also update the column so that subsequent ambiguous bases can be found
							# as these are indexed by column.index(j) so the same entry is updated
							# repeatedly instead of the next one being updated
							column = column[:column.index(j)] + BB + column[column.index(j) + 1:]
		if verbose:
			if monobaseColumns == 0:
				out = 0
			else:
				out = monobaseCount / float(monobaseColumns)
			return ['Ambiguous bases corrected: %i' % (monobaseCount),
				'Columns containing ambiguous bases: %i (%3.2f%%)' % (monobaseColumns, monobaseColumns/float(self.seqLength()[0][1])*100),
				'Ambiguous bases per corrected column: %3.2f' % (out)]
	# ----------------

	def blunt(self, left=True, right=True):
		'''Remove columns from multiple sequence alignment such that no sequence starts or ends with a gap'''
		# method to trim off tatty ends (make blunt ends)
		# required for phylogenetic analysis
		# but not required when joining fragments together, as more bases should be included in this case
		## import operator
		## X = sorted(self.seq, key=operator.itemgetter('seq'))	# sort according to actual sequence data; slow for large sequences? memory requirements?

		if left:
			longest = 0
			for i in range(len(self.seq)):	# process all sequences
				for j in range(longest, self.seqLength()[0][1]):
					if self.seq[i]['seq'][j] != GAP:
						longest = j
						break

			leftCount = longest

			for i in range(len(self.seq)):
				self.seq[i]['seq'] = self.seq[i]['seq'][longest:]

		else:
			leftCount = None

		if right:
			longest = 0
			for i in range(len(self.seq)):	# process all sequences
				for j in range(self.seqLength()[0][1] - longest - 1, 0, -1):
					if self.seq[i]['seq'][j] != GAP:
						longest = self.seqLength()[0][1] - j - 1
						break

			rightCount = longest

			longest = self.seqLength()[0][1] - longest
			for i in range(len(self.seq)):
				self.seq[i]['seq'] = self.seq[i]['seq'][:longest]

		else:
			rightCount = None

		return [leftCount, rightCount]
	# ---------------

	def seqSerotype(self, startMotif='ATGGAGAACAT'):
		'''Return serotype according to Purdy 2007'''
		# Requires aligned sequences
		# Using list.index() to extract each 'find' result can be used to get around this later if necessary

		import Bio.SeqUtils

		StartS = self.find(startMotif)[0][2]

		seroMotifs = self.extract('122,160,127,159,140', aain=True, aaout=True,  mapping = -StartS + 2)
		seroNucs   = self.extract('122,160,127,159,140', aain=True, aaout=False, mapping = -StartS + 2)
		out = []

		for i in seroMotifs:
			serotyped = False
			# print i[0][:20], i[1],
			for j in SEROPATTERNS:
				found = re.search(j, i[1])
				if found:
					# print j, SEROPATTERNS[j]
					serotype = SEROPATTERNS[j]
					serotyped = True
			if not serotyped:
				if i[1][0] in SEROGROUP:
					# print j[0], SEROGROUP[j[0]]	# at least "ad" or "ay"  --  i[1][0] is the first amino acid in the pattern
					serotype = SEROGROUP[i[1][0]]
				else:
					# print "Unknown"
					serotype = "Unknown"

			three = ''
			for j in i[1]:	# each letter in the motif
				three += Bio.SeqUtils.seq3(j)

			out.append([i[0], serotype, i[1], three, seroNucs[seroMotifs.index(i)][1]])
		return out
	# ---------------

	def mutationFinder(self, reference, startMotif, stopMotif, gene, aa=True, nucleotideMapping=0):
		'''Return mutations found in self using given object reference as reference sequence'''
		# patterns are regular expressions
		# aa=True means return amino acid positions
		# currently only the first sequence in the reference object is processed
		self.seqCase()
		reference.seqCase()

		prefixList = {'P' : 'rt', 'S': 's', 'BCP' : '', 'X': 'x', 'C' : 'c'}

		if gene not in prefixList:
			prefixList[gene] = ''	# ensure lookup later works

		startReference = re.search(startMotif, reference.seq[0]['seq'])
		stopReference  = re.search(stopMotif, reference.seq[0]['seq'])

		if startReference == None or stopReference == None:
			return None			# if either motif is not found in the reference sequence, abort

		startReference = startReference.start()
		stopReference  = stopReference.start()

		out = []	# will store the list out found mutations: ['ID1', ['A123C', 'C456G']]

		for theSeq in self.seq:
			startQuery = re.search(startMotif, theSeq['seq'])
			stopQuery  = re.search(stopMotif, theSeq['seq'])

			if startQuery == None or stopQuery == None:
				out.append([theSeq['id'], []])		# return (iterable) empty list rather than (non-iterable) "None"
				continue		# start or stop motifs are not found; abort this sequence

			startQuery = startQuery.start()
			stopQuery  = stopQuery.start()

			mapping = startReference - startQuery

			outQuery = ''
			outReference = ''
			mutationList = []	# temporarily store the list of mutations

			for i in range(startQuery, stopQuery):
				outQuery += theSeq['seq'][i]
				outReference += reference.seq[0]['seq'][i + mapping]

			# should avoid building and translating the reference sequence each time ...

			if aa:
				#import Bio.Seq

				#transReference = Bio.Seq.translate(outReference)
				#transQuery = Bio.Seq.translate(outQuery)

				#for i in range(len(transQuery)):
				#	if transQuery[i] != transReference[i]:
				#		mutationList.append( '%s%s%i%s' % (prefixList[gene], transReference[i], i+1, transQuery[i]) )	# startReference + (i * 3) = nucleotide position

				import Bio.Seq

				for i in range(0, len(outQuery), 3):
					rr = outReference[i:i+3]
					qq = outQuery[i:i+3]
					if rr != qq:
						mutationList.append((i+startQuery, rr, qq, Bio.Seq.translate(rr), Bio.Seq.translate(qq), i/3+1))	# position output is the actual position in the sequence, not the amino acid position



			else:
				for i in range(len(outReference)):
					if outQuery[i] != outReference[i]:
						mutationList.append('%s%s%i%s' % (prefixList[gene], outReference[i], i+nucleotideMapping, outQuery[i]) )	# 1741

			out.append([theSeq['id'], mutationList])


		return out
	# ---------------

	def seqNames(self):
		out = []
		for i in self.seq:
			out.append(i['id'])
		return out
	# ---------------

	def seqGenotype(self, s_start=S_START, b_start=B_START):
		'''Return serotype according to Kramvis 2008'''
		# adapted from stand-alone 'genotyper.py'

		import Bio.SeqUtils

		out = []
		for i in self.seq:         # iterate over each sequence searching each one; not using S.find which searches all
			ss = re.search(S_START, i['seq'])
			bb = re.search(B_START, i['seq'])

			if ss != None:
				ss = ss.start()+1
			if bb != None:
				bb = bb.start()+1



			if ss != None:
				seroMotif = self.extract('122,160,127,159,140', aain=True, aaout=True,  mapping = -ss + 2, inData=[i])[0]
				seroNucs   = self.extract('122,160,127,159,140', aain=True, aaout=False, mapping = -ss + 2, inData=[i])[0]

				serotyped = False
				for k in SEROPATTERNS:
					found = re.search(k, seroMotif[1])
					if found:
						serotype = SEROPATTERNS[k]
						serotyped = True
				if not serotyped:
					if seroMotif[1][0] in SEROGROUP:
						serotype = SEROGROUP[seroMotif[1][0]]
					else:
						serotype = '??'

				three = ''
				for k in seroMotif[1]:	# each letter in the motif
					three += Bio.SeqUtils.seq3(k)

				serotype = ([seroMotif[0], serotype, seroMotif[1], three, seroNucs[1]])
			else:
				serotype = '??'

			if bb != None:
				bcpMotif = self.extract('1802,1803,1809-1812,1817-1819,1858,1888', aain=False, aaout=False, mapping = 1814-bb+1, inData=[i])[0]
			else:
				bcpMotif = '??'

			if bb != None and ss != None:

				final = serotype[1] + bcpMotif[1]

				found = False
				for j in GENOTYPES:
					found = re.search(j, final)

					if found:
						genotype = (i['id'][:self.truncate], GENOTYPES[j], final, serotype[1], serotype[2], serotype[3], serotype[4], len(i['seq']))
						break

				if not found:
					genotype = (i['id'][:self.truncate], '??', final, serotype[1], serotype[2], serotype[3], serotype[4], len(i['seq']))

			else:
				genotype = (i['id'][:self.truncate], '??', '??', '??' ,'??', '??', '??', len(i['seq']))

			out.append(genotype)

		return (out)

	# ---------------

	def DEPRECATED_groupDistribution_DEPRECATED(self, grouping, threshold=0.05, aa=False):
		'''Return chi-squared ViPR/BRC p-value for all columns using provided grouping'''
		# Passing parameter for all columns is better, as then the splitting into groups is done only once,
		# rather than repeating this for each column. Typical usage would be to process all columns in a sequence.

		# Differences from the ViPR tool: this one allows for more groups and *requires* aligned sequences
		# how to identify a sequence is protein not nucleotide; conceivable that a *protein* column
		# contain only A, C, G and T amino acids...?

		# Web interface: file to upload -- threshold -- amino acid -- grouping list for 12 groups
		# Sequences *must* be aligned -- this is better, because the sequence data can be checked/curated before
		# Otherwise, if an alignment is done and then the data is scanned for significant columns, these results
		# would depend on the (unchecked) alignment

		# MRT is used to examine distribution at loci of interest which are known before
		# This tool can focus on *which* positions are of interest and then these can be given to the MRT to report on.
		# The two tools are used together -- this one then the MRT

		# Must abort if groups < 2
		# Stop checking if no variation in a column (rather than ending up with degrees of freedom = 0)

		# Algorithm for dynamic totals for amino acids:
		# - Start as for nucleotides, totalling each into it's own column
		# - Afterwards, total the number of non-zero columns

		import scipy, scipy.stats, numpy

		emptyColumn = 0.01			# used instead of zero in empty columns

		groups = len(grouping)
		gg = []
		for i in grouping:			# each group is expanded into list gg
			gg.append(parseList(i))		# ['1-10', '11,13'] becomes [[[1,10]], [[11,11],[13,13]]]

		g = []					# which then becomes [(1,2,3,4,5,6,7,8,9,10), (11,13)]   (1-indexed here)
		for i in gg:				# stores the grouping as a list of tuples (0-indexed)
			tt = ()
			for j in i:
				for k in range(j[0], j[1]+1):
					tt += (k-1,)	# k-1: convert from 1-indexed (provided grouping) to 0-indexed (required)
			g.append(tt)

		residueCols = {}	# stores the position (column number) for each residue; A in column 0, B in column 1, etc.
		cc = 0

		if not aa:
			residueList = BASES
			residueCols[GAP] = len(BASES)	# gaps are stored in the last (right-most) column
		else:
			residueList = AMINOACIDS
			residueCols[GAP] = len(AMINOACIDS)

		for i in residueList:
			residueCols[i] = cc
			cc += 1

		for i in range(len(self.seq[0]['seq'])):		# each column in the sequence
			if not aa:
				residueCount = len(BASES) + 1		# for nucleotides; dynamically determined for amino acids below

			observed = numpy.zeros((groups, len(residueCols)), dtype=float)
			expected = numpy.zeros((groups, len(residueCols)), dtype=float)
			gapsPresent = False		# flag; if at least one gap, increase degrees of freedom
			for j in g:			# each group
				for k in j:		# each sequence in each group
					if self.seq[k]['seq'][i] in residueList + [GAP]:	# ignore non-bases or non-gaps
						observed[g.index(j), residueCols[self.seq[k]['seq'][i]]] += 1
					if self.seq[k]['seq'][i] == GAP:
						gapsPresent = True
					# looks up each residue in the dictionary to determine the column number
					# row number (group number) determined as follows:
					# j could hold (1,2,3) for group '1' for example, so index returns the group number

			if not aa:
				degreesFreedom = (groups - 1) * (4 - 1)		# nt: dof fixed at 4
			else:
				residueCount = 0
				for j in range(len(AMINOACIDS)+1):
					if observed.sum(axis=0)[j] > 0:		# determine which residues have a non-zero total
						residueCount += 1
				degreesFreedom = (groups - 1) * (residueCount - 1)

			if not gapsPresent:

				if not aa:
					residueCount -= 1
	     				observed = numpy.delete(observed, 4, 1)		# delete column (axis 1) number 4
					expected = numpy.delete(expected, 4, 1)		# delete column (axis 1) number 4
					# for amino acids, the gap is just another residue which is included in the preceeding code
			else:
				degreesFreedom =+ 1	# if gaps present, add one to degrees of freedom

			for j in range(groups):
				for k in range(residueCount):
					if observed[j, k] == 0:
						observed[j, k] = emptyColumn
			for j in range(groups):
				expected[j] = (observed.sum(axis=0) * observed[j].sum()) / observed.sum()

			chisq = ((observed - expected) ** 2) / expected

			if i == 27:
				print "observed", observed.round(decimals=3)
				print "expected", expected.round(decimals=3)
				print residueCount, "residue count"
				print degreesFreedom, "degrees freedom"
				print "CALCULATION", ((observed - expected) ** 2) / expected
				print "CHI SQUARED", chisq.sum()

			CHI = scipy.stats.chisqprob(chisq.sum(), degreesFreedom)
			if CHI <= threshold:
				print 'position %04i, p=%3.8f, chi-squared=%3.2f, DOF=%i' % (i+1, CHI, chisq.sum(), degreesFreedom)

	# ---------------

	def parseReGroups(self, grouping):
		import re
		listID = {}
		c = 0
		for i in self.seq:
			listID[c] = i['id']                             # ID stored as dictionary value rather than key; the key is the position
			c += 1

		numSeqs = len(listID)
		c = 0
		groups = []
		for i in grouping:
			thisGroup = ''
			p = 0
			while (len(listID) > 0) and (p <= numSeqs):            # while there are items left in the dictionary and we haven't gone further than the initial number of sequences
				if p in listID:
					match = re.search(i, listID[p])         # search for the group ("D4") in each of the keys of the dictionary
					if match:
						c += 1
						thisGroup = thisGroup + str(p+1) + ','		# positions are 1-indexed
						del listID[p]
					p += 1
				else:
					p += 1
			groups.append(thisGroup[:-1])                           # strip the trailing comma and add to the groups list; tuple
		# ----- i -----
		return groups
	# ----- parseReGroups -----

	def groupDistribution(self, grouping, threshold=0.05, aa=False, reGroups=False):

		# Assumptions: input sequences are aligned and of identical length
		# Non-nucleotide and/or non-amino acid bases
		# Completely ignore completely conserved columns

		import scipy, scipy.stats, numpy

		emptyColumn = 0.01			# used instead of zero in empty columns

		if reGroups:
			grouping = self.parseReGroups(grouping)		# convert to standard grouping format: ['1,2-10,12', '11,15-19']

		groups = len(grouping)
		gg = []
		for i in grouping:			# each group is expanded into list gg
			gg.append(parseList(i))		# ['1-10', '11,13'] becomes [[[1,10]], [[11,11],[13,13]]]
		g = []					# which then becomes [(1,2,3,4,5,6,7,8,9,10), (11,13)]   (1-indexed here)
		for i in gg:				# stores the grouping as a list of tuples (0-indexed)
			tt = ()
			for j in i:
				for k in range(j[0], j[1]+1):
					tt += (k-1,)	# k-1: convert from 1-indexed (provided grouping) to 0-indexed (required)
			g.append(tt)			# g is a list of tuples of each group

		if aa:
			residueList = AMINOACIDS[:]	# copy list
		else:
			residueList = BASES[:]		# copy list
		residueList += [GAP]			# add the GAP so that we can check to see if any are present later

		for i in range(len(self.seq[0]['seq'])):			# each column
		# for i in range(1):
			residues = []						# list (for column) holding dictionary (for each group)
			for j in g:						# each group
				rr = {}						# make fresh each time; 'copying' didn't seem to work
				for rrLoop in residueList:
					rr[rrLoop] = 0

				for k in j:					# each sequence (in each group)
					theResidue = self.seq[k]['seq'][i]
					if theResidue not in residueList:
						theResidue = GAP
					rr[theResidue] += 1			# increment each residue
				# ----- end sequence -----
				residues.append(rr)				# append the dictionary (rr) for the group to the list (residues) for the column
			# ----- end group -----

			# Check here if column is completely conserved?

			if aa:
				# all amino acids, or a gap, are considered as a residue in the column; total them up into reisdueCount
				rr = {}						# make fresh each time; 'copying' didn't seemd to work
				for rrLoop in residueList:
					rr[rrLoop] = False

				residueCount = 0
				residuesFound = []
				for j in residues:
					for k in j:
						if (j[k] > 0) and (not rr[k]):	# if there is >0 residue occurences and we've not seen this residue before...
							rr[k] = True		# ...then mark it as seen...
							residueCount += 1	# ...and add to the residueCount
							residuesFound.append(k)
			# ----- aa -----
			else:
				# all four nucleotides are always automatically included; if at least one gap is present, then it is also included
				residueCount = 4
				residuesFound = BASES[:]

				for j in residues:
					if j[GAP] > 0:
						residueCount = 5
						residuesFound = residueList	# includes GAP
						break

				for j in residues:
					for k in j:
						if j[k] == 0:
							j[k] = emptyColumn	# all four nucleotides are included; those which are absent are scored at 'emptyColumn'
			# ----- not aa -----

			# print 'residues', residues
			# print 'residuecount', residueCount
			# print 'residuesFound', residuesFound

			observed = numpy.zeros((groups, residueCount), dtype=float)
			expected = numpy.zeros((groups, residueCount), dtype=float)

			degreesFreedom = (groups - 1) * (residueCount - 1)

			residuesFound.sort()

			for j in range(groups):
				for k in range(residueCount):
					observed[j, k] = residues[j][residuesFound[k]]

			for j in range(groups):
				expected[j] = (observed.sum(axis=0) * observed[j].sum()) / observed.sum()

			chisq = ((observed - expected) ** 2) / expected

			CHI = scipy.stats.chisqprob(chisq.sum(), degreesFreedom)
			if CHI <= threshold:
				print 'position %04i, p=%3.8f, chi-squared=%3.2f, DOF=%i' % (i+1, CHI, chisq.sum(), degreesFreedom)

			# print "observed", observed.round(decimals=3)
			# print "expected", expected.round(decimals=3)
			# print residueCount, "residue count"
			# print degreesFreedom, "degrees freedom"
			# print "CALCULATION", ((observed - expected) ** 2) / expected
			# print "CHI SQUARED", chisq.sum()

		# ----- end column -----

	# ---------------

	def wt2x2(self, grouping, threshold=0.1, reGroups=False):
		'''Return position, wild-type, Fisher's and Chi for significant positions in 2x2 contingency table'''

		import math, scipy.stats

		if reGroups:
			grouping = climb.parseReGroups(grouping)		# convert to standard grouping format: ['1,2-10,12', '11,15-19']

		groups = len(grouping)
		gg = []
		for i in grouping:			# each group is expanded into list gg
			gg.append(parseList(i))		# ['1-10', '11,13'] becomes [[[1,10]], [[11,11],[13,13]]]
		g = []					# which then becomes [(1,2,3,4,5,6,7,8,9,10), (11,13)]   (1-indexed here)
		for i in gg:				# stores the grouping as a list of tuples (0-indexed)
			tt = ()
			for j in i:
				for k in range(j[0], j[1]+1):
					tt += (k-1,)	# k-1: convert from 1-indexed (provided grouping) to 0-indexed (required)
			g.append(tt)			# g is a list of tuples of each group

		out = []

		length = len(self.seq[0]['seq'])
		for i in range(0, length):              # each position in the sequence
			d = self.baseDistribution(str(i+1), set1=AMINOACIDS, set2=[GAP])
			# print d
			high = 0
			highRes = GAP
			for j in d[0]:
				if int(d[0][j]) > high:
					high = int(d[0][j])
					highRes = j

			# iterate over g[0] for group 1; iterate over g[1] for group 2
			mutantList1 = ''
			mutantList2 = ''

			tempC = 0
			for j in g[0]:
				if self.seq[j]['seq'][i] == highRes:
					tempC += 1
				else:
					mutantList1 += self.seq[j]['seq'][i]
			group1 = (tempC, len(g[0]) - tempC)

			tempC = 0
			for j in g[1]:
				if self.seq[j]['seq'][i] == highRes:
					tempC += 1
				else:
					mutantList2 += self.seq[j]['seq'][i]
			group2 = (tempC, len(g[1]) - tempC)

			# print group1, group2
			# 2x2 contingency "shortcut" variables
			a = group1[0]
			b = group1[1]
			c = group2[0]
			d = group2[1]

			FET = (math.factorial(a+b) * math.factorial(c+d) * math.factorial(a+c) * math.factorial(b+d)) / float(math.factorial(a) * math.factorial(b) * math.factorial(c) * math.factorial(d) * math.factorial(a+b+c+d))

			# if FET < 0.1:
			# 	print "position %03i with wild-type %s has p=%4.3f (Fisher's Exact)" % (i+1, highRes, FET)
			# 	if verbose:
			# 		print "%05s %05s %05s %05s" % ("grp", "wt", "mut", "tot")
			# 		print "%05s %5i %5i %5i" % ("1", group1[0], group1[1], group1[0] + group1[1])
			# 		print "%05s %5i %5i %5i" % ("2", group2[0], group2[1], group2[0] + group2[1])
			# 		print "%05s %5i %5i %5i" % ("tot", group1[0] + group2[0], group1[1] + group2[1], group1[0] + group1[1] + group2[0] + group2[1])
			# 		print

			if (a > 5) and (b > 5) and (c > 5) and (d > 5):
				CHI = (((a*d - b*c) ** 2) * (a+b+c+d)) / float( (a+b) * (c+d) * (b+d) * (a+c))
				CHIp = scipy.stats.chisqprob(CHI, 1)    # 1 dof
				# if CHIp < 0.1:
				# 	print "position %03i with wild-type %s has p=%4.3f (Chi Square = %4.3f)" % (i+1, highRes, CHIp, CHI)
			else:
				CHIp = 99	# prevent including this in the returned output below

			# returns position, wild-type, FET, CHIp, group1, group2
			if (FET <= threshold) or (CHIp <= threshold):
				out.append([i+1, highRes, FET, CHIp, group1, group2, ''.join(sorted(mutantList1)), ''.join(sorted(mutantList2))])

		return out
	# ---------------

	def gvt(self, grouping, reGroups=False):
		'''Return variation at each locus for two groups'''

		# Effectively, a base distribution method for two groups

		if reGroups:
			grouping = climb.parseReGroups(grouping)		# convert to standard grouping format: ['1,2-10,12', '11,15-19']

		groups = len(grouping)
		gg = []
		for i in grouping:			# each group is expanded into list gg
			gg.append(parseList(i))		# ['1-10', '11,13'] becomes [[[1,10]], [[11,11],[13,13]]]
		g = []					# which then becomes [(1,2,3,4,5,6,7,8,9,10), (11,13)]   (1-indexed here)
		for i in gg:				# stores the grouping as a list of tuples (0-indexed)
			tt = ()
			for j in i:
				for k in range(j[0], j[1]+1):
					tt += (k-1,)	# k-1: convert from 1-indexed (provided grouping) to 0-indexed (required)
			g.append(tt)			# g is a list of tuples of each group

		out = []

		length = len(self.seq[0]['seq'])
		for i in range(0, length):              # each position in the sequence

			mutantList1 = {}	# dictionaries to store the number of occurrences of each reside in each group
			mutantList2 = {}

			# iterate over g[0] for group 1; iterate over g[1] for group 2
			for j in AMINOACIDS + ['U', 'O', 'B', 'Z' ,'J', '*']:
				mutantList1[j] = 0
				mutantList2[j] = 0
			for j in BASES + [GAP] + ['N', '?', 'X']:
				mutantList1[j] = 0
				mutantList2[j] = 0

			for j in g[0]:
				mutantList1[self.seq[j]['seq'][i]] += 1

			for j in g[1]:
				mutantList2[self.seq[j]['seq'][i]] += 1

			out.append([mutantList1, mutantList2])

		return (out, g)		# return the list of distributions and the grouping
	# ---------------

	def deepThreshold(self, prErr=0.01, power=0.05):
		import scipy.stats
		df = 1
		def chi2(o, e):
			return ((o-e)**2)/(e)
		# --------------
		chiThreshold = scipy.stats.chi2.isf(power, df)
		Expected = int(round(prErr * self.seqCount()))	# float(prErr * self.seqCount())
		if Expected == 0:
			Expected = 1
		Observed = Expected+1
		while chi2(Observed, Expected) < chiThreshold:
			Observed += 1
		return (Expected, Observed)
	# --------------

	def deepColumns(self, countThreshold, numColumns=0):
		# count residues in column which are not majority
		# any above threshold mean that the position is interesting
		errorsPerRow = []
		for i in range(self.seqCount()):
			errorsPerRow.append(0)
		c = 1
		interestingColumns = []
		if numColumns == 0:
			numColumns = len(self.seq[0]['seq'])
		for i in self.baseDistribution('1-%s' % (str(numColumns))):
			maxCount = 0
			maxResidue = ''
			ignoreBases = []
			for j in (BASES + NONBASES):
				if i[j] > maxCount:
					maxCount = i[j]
					maxResidue = j
			for j in (BASES + NONBASES):
				if (j <> maxResidue) and (i[j] >= countThreshold):
					interestingColumns.append(c)
					break
			for j in (BASES + NONBASES):
				if i[j] < countThreshold:
					ignoreBases.append(j)

			r = 0
			for j in self.seq:
				if j['seq'][c-1] in ignoreBases:
					errorsPerRow[r] += 1
				r += 1

			c += 1
		return (interestingColumns, errorsPerRow)
	# --------------

        def SRID(self, searchString, replaceString):
                import re
                import copy
                # >gi|87295395|gb|DQ315783.1| Hepatitis B virus isolate EI03188, complete genome
                # re.sub(r'.+[bj]\|(?P<n>.+)\.1\|.+', r'\g<n>', x)
                Temp = copy.deepcopy(self)
                for i in Temp.seq:
                        i['id'] = re.sub(searchString, replaceString, i['id'])
                return Temp
                # S.writeFASTA(sys.argv[2])
	# --------------

        def SplitID(self, splitList):
                import re
                import copy
                # splitList = ['Y18855','Y18856','Y18857','Y18858','DQ315778','FJ023675','FJ023659']
                Temp1 = copy.deepcopy(self)
                Temp2 = copy.deepcopy(self)

                Temp1.seqRemoveByID(splitList, retain=True, verbose=False, escape=True)
                Temp2.seqRemoveByID(splitList, retain=False, verbose=False, escape=True)

                return (Temp1, Temp2)
	# --------------

        def SplitLength(self, expression):
                import copy
                Temp1 = copy.deepcopy(self)
                Temp2 = copy.deepcopy(self)
                matchList = []
                for i in self.seq:
                        l = len(i['seq'])
                        if eval(expression):
                                matchList.append(i['id'])

                Temp1.seqRemoveByID(matchList, retain=True, verbose=False, escape=True)
                Temp2.seqRemoveByID(matchList, retain=False, verbose=False, escape=True)

                return (Temp1, Temp2)
	# --------------

        def SRseq(self, searchString, replaceString):
                import copy
                Temp = copy.deepcopy(self)
                Temp.seqSR(searchString, replaceString)
                return Temp
        # --------------

        def SplitMotif(self, queryLoci, queryMotif):
                import copy
                Temp1 = copy.deepcopy(self)
                Temp2 = copy.deepcopy(self)
                matchList = []
                t = self.extract(queryLoci)
                for i in range(0, self.seqCount()):
                        if re.match(queryMotif, t[i][1]) is not None:
                                matchList.append(self.seq[i]['id'])

                Temp1.seqRemoveByID(matchList, retain=True, verbose=False, escape=True)
                Temp2.seqRemoveByID(matchList, retain=False, verbose=False, escape=True)

                return (Temp1, Temp2)
        # --------------

if __name__ == '__main__':
	print motd()

