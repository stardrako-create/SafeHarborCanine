from review_sv_pair_evidence_20260922 import classify,CHROM,BOUNDS
import copy
l=dict(mate=1,chrom=CHROM,start0=BOUNDS[0],mapq=60,reverse=False,paired=True,unmapped=False,secondary=False,supplementary=False,qcfail=False,duplicate=False,mate_unmapped=False,mate_chrom=CHROM,mate_start0=BOUNDS[1],mate_reverse=True)
r=dict(l,mate=2,start0=BOUNDS[1],reverse=True,mate_start0=BOUNDS[0],mate_reverse=False)
assert classify(l,r)=='both_MAPQ30_inward'
assert classify(l,dict(r,mapq=29))=='low_MAPQ'
assert classify(dict(l,mapq=29),r)=='low_MAPQ'
assert classify(l,dict(r,duplicate=True))=='ineligible_record'
assert classify(l,None)=='missing_or_ambiguous_mate_record'
assert classify(l,dict(r,mate_start0=1))=='nonreciprocal_mate_fields'
assert classify(l,dict(r,mate=1))=='mate_identity_mismatch'
assert classify(dict(l,reverse=True,mate_reverse=False),dict(r,reverse=False,mate_reverse=True))=='both_MAPQ30_other_orientation'
assert classify(l,dict(r,mate_reverse=True))=='inconsistent_mate_strand_fields'
print('9 targeted assertions passed')
