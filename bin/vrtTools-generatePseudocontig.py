#!/usr/bin/env python

import sys
import subprocess as sp

if len(sys.argv)!=3:
	print "::usage: %s <fasta_reference> <fasta_draft>" % sys.argv[0]
	sys.exit()

ref=sys.argv[1]
draft=sys.argv[2]

cmd="python ~/sw/contiguator/2.7/CONTIGuator.py -r %s -c %s" % (ref,draft)
sp.call(cmd,shell=True)
