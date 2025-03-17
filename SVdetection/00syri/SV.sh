minimap2 -x asm5 -t 20 -c --eqx -o FDA.FDB.paf $FDA $FDB 
syri -c FDA.FDB.filter.paf -r $FDA -q $FDB -F P -f --nc 15 --cigar
