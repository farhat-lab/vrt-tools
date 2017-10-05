#!/usr/bin/env python

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

def read_snps_generate_vrt(mummer_out,ref_seq,pcontig_seq,file_out):
	with open(file_out,"w") as outf:
		with open(mummer_out,"r") as f:
			for i in xrange(5):
				next(f)

			case1_ref=""
			case1_alt=""
			case1_pos=""
			case1_strand=""

			case2_ref=""
			case2_alt=""
			case2_pos=""
			case2_strand=""
			counter=0
			for line in f:
				counter=counter+1
				print "\r--analyzing line %d\t\t\t" % counter,
				mum_data=re.split("\s+",line)
				#outf.write("\t".join(mum_data)+"\n")
				pos_ref=int(mum_data[1])
				pos_alt=int(mum_data[4])
				ref=mum_data[2]
				alt=mum_data[3]
				if((ref!=".") and (alt!=".")):
					if(case1_pos!=""):
						if(case1_strand=="-"):
							bseq=Seq(case1_alt)
							outf.write("INS\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case1_pos),case1_ref,str(bseq.complement()),len(str(bseq.complement()))-len(case1_ref))) ##complement!
						else:
							outf.write("INS\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case1_pos),case1_ref,case1_alt,len(case1_alt)-len(case1_ref))) ##complement!
						case1_ref=""
						case1_alt=""
						case1_pos=""
						case1_strand=""
						outf.write("SNP\t%d\t%s\t%s\t1\tNA\tNA\n" % (int(pos_ref),ref,alt))
					elif(case2_pos!=""):
						if(case2_strand=="-"):
							bseq=Seq(case2_alt)
							outf.write("DEL\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case2_pos),case2_ref,str(bseq.complement()),len(case2_ref)-len(str(bseq.complement())))) ##complement!
						else:
							outf.write("DEL\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case2_pos),case2_ref,case2_alt,len(case2_ref)-len(case2_alt)))
						case2_ref=""
						case2_alt=""
						case2_pos=""
						case2_strand=""
						outf.write("SNP\t%d\t%s\t%s\t1\tNA\tNA\n" % (int(pos_ref),ref,alt))
					else:
						outf.write("SNP\t%d\t%s\t%s\t1\tNA\tNA\n" % (int(pos_ref),ref,alt))
				if((ref==".") and (case1_pos=="") and (alt!=".")):
					if(case2_pos!=""):
						if(case2_strand=="-"):
							bseq=Seq(case2_alt)
							outf.write("DEL\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case2_pos),case2_ref,str(bseq.complement()),len(case2_ref)-len(str(bseq.complement())))) ##complement!
						else:
							outf.write("DEL\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case2_pos),case2_ref,case2_alt,len(case2_ref)-len(case2_alt)))
					
						case2_ref=""
						case2_alt=""
						case2_pos=""
						case2_strand=""
					case1_pos=(pos_ref)

					if(mum_data[16]=="-1"):

						case1_ref=list(ref_seq)[(pos_ref-1)]
						case1_alt=list(pcontig_seq)[(pos_alt)]
						case1_strand="-"
					else:
						case1_ref=list(ref_seq)[(pos_ref-1)]
						case1_alt=list(pcontig_seq)[(pos_alt-2)]
				if((ref==".") and (case1_pos!="") and (alt!=".")):
					case1_alt=case1_alt+alt
				if((ref!=".") and (case2_pos=="") and (alt==".")):
					if(case1_pos!=""):
						if(case1_strand=="-"):
							bseq=Seq(case1_alt)
							outf.write("INS\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case1_pos),case1_ref,str(bseq.complement()),len(str(bseq.complement()))-len(case1_ref))) ##complement!
						else:
							outf.write("INS\t%d\t%s\t%s\t%d\tNA\tNA\n" % (int(case1_pos),case1_ref,case1_alt,len(case1_alt)-len(case1_ref))) ##complement!
						case1_ref=""
						case1_alt=""
						case1_pos=""
						case1_strand=""
					case2_pos=(pos_ref-1)

					if(mum_data[16]=="-1"):
						case2_ref=list(ref_seq)[(pos_ref-2)]
						case2_alt=list(pcontig_seq)[(pos_alt)]
						case2_strand="-"
					else:
						case2_ref=list(ref_seq)[(pos_ref-2)]
						case2_alt=list(pcontig_seq)[(pos_alt-1)]
				if((ref!=".") and (case2_pos!="") and (alt==".")):
					case2_ref=case2_ref+ref


def main():
	if len(sys.argv)!=5:
		print "::usage: %s <snps_file> <fasta_ref> <fasta_pseudocontig> <prefix_out_files>" % sys.argv[0]
		sys.exit()

	mummer_out=sys.argv[1]
	fasta_ref=sys.argv[2]
	fasta_pcontig=sys.argv[3]
	file_out=sys.argv[4]

	ref = SeqIO.read(fasta_ref, "fasta")
	p_contig= SeqIO.read(fasta_pcontig, "fasta")

	read_snps_generate_vrt(mummer_out,str(ref.seq),str(p_contig.seq),file_out)
	print "\nBye!\n"


if __name__ == "__main__":
	main()



