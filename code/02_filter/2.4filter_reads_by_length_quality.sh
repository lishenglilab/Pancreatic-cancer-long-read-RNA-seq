# To filter reads that are EITHER <200bp OR have Q<7 

#!/bin/bash


cd /public/nanopore_cDNA/fastq/

find . -type f -name "filtered_*.fastq" -print0 | \
parallel -0 -j 20 '
    file_dir="{//}"       
    file_name="{/}"       
    base_name="{/.}"      

    cd "$file_dir" && \
    NanoFilt -l 200 -q 7 "$file_name" > "${base_name}_nanofilt.fastq"
'
