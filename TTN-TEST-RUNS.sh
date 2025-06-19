unzip TTN-TEST-RUNS.zip
#get ODB database for BUSCO
wget https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz
tar xvf eukaryota_odb10.2024-01-08.tar.gz
#TEST DATA TITIN
#run quick test with different evidence (prot, prot+rna) and with or without BUSCO
#switch off EGGNOG for speed (-E 1)
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+r-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa | bash > RUN-TEST1-p+r-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa | bash > RUN-TEST1-p-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+r+b-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa -b eukaryota_odb10 | bash > RUN-TEST1-p+r+b-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+b-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -b eukaryota_odb10 | bash > RUN-TEST1-p+b-noegg.log 2>&1
#run quick test with different evidence (prot+gtf, prot+rna+gtf) and with or without BUSCO
#switch off EGGNOG for speed (-E 1)
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+r+g-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa -g TTN_canary.gtf | bash > RUN-TEST1-p+r+g-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+g-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -g TTN_canary.gtf | bash > RUN-TEST1-p+g-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+r+b+g-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa -b eukaryota_odb10 -g TTN_canary.gtf | bash > RUN-TEST1-p+r+b+g-noegg.log 2>&1
./scripts/HANNO.v0.6.pl -E 1 -t 8 -d RUN-TEST1-p+b+g-noegg -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -b eukaryota_odb10 -g TTN_canary.gtf | bash > RUN-TEST1-p+b+g-noegg.log 2>&1
#run quick test with different evidence (prot, prot+rna) and with or without BUSCO
./scripts/HANNO.v0.6.pl  -t 8 -d RUN-TEST1-p+r   -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa | bash > RUN-TEST1-p+r.log 2>&1
./scripts/HANNO.v0.6.pl  -t 8 -d RUN-TEST1-p     -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa | bash > RUN-TEST1-p.log 2>&1
./scripts/HANNO.v0.6.pl  -t 8 -d RUN-TEST1-p+r+b -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -r TTN-mRNA-zebrafinch.fa -b eukaryota_odb10 | bash > RUN-TEST1-p+r+b.log 2>&1
./scripts/HANNO.v0.6.pl  -t 8 -d RUN-TEST1-p+b   -a TTN-genomic-canary.fa -p TTN-AA-zebrafinch.fa -b eukaryota_odb10 | bash > RUN-TEST1-p+b.log 2>&1
