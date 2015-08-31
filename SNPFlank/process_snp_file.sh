#! /bin/bash

# Remove bad groups
grep -v -f bad_group_patterns ckqny.Groups.txt >ckqny.Groups.filtered.txt

# Refactor to bed format
awk 'BEGIN{OFS="\t";FS="\t"}; $1~/^group/{print $5,$6,$6,$2,$7,"+"};' ckqny.Groups.filtered.txt >ckqny.Groups.bed

# Get rsIDs
cut -f4 ckqny.Groups.bed | grep "^rs" > all_SZ_snps_rs.txt
