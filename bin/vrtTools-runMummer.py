#!/usr/bin/env python

import sys
import subprocess as sp

if len(sys.argv)!=4:
	print "::usage: %s <fasta_reference> <fasta_pseudocontig> <prefix_out_files>" % sys.argv[0]
	sys.exit()

ref=sys.argv[1]
p_contig=sys.argv[2]
prefix=sys.argv[3]

cmd="nucmer --prefix=%s %s %s" % (prefix,ref,p_contig)
sp.call(cmd,shell=True)
cmd="show-snps -Clr %(p)s.delta > %(p)s.snps" % {"p":prefix}
sp.call(cmd,shell=True)
