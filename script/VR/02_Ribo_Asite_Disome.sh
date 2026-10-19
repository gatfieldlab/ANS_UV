### LICENCE
# Copyright (c) 2026 Virginie Ricci

# This file is part of ANS_UV (https://github.com/gatfieldlab/ANS_UV/).

# ANS_UV includes free scripts: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.

# ANS_UV is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this script. 
# If not, see <https://www.gnu.org/licenses/>.
### LICENCE 

module load chip-seq/1.5.5 # https://epd.expasy.org/chipseq/
module load samtools/1.19.2 

SHIFT=45

# Asite offset: + 45 for disome
for i in /path_to/ANS_HD_UV_LD/RSEM/sub_bam/PD*.filtered.bam; do
        filename=$(basename $i .transcript.sorted.filtered.bam)
        echo $filename

        samtools view $i | awk -v "sh=$SHIFT" 'BEGIN{FS=OFS="\t"}{id=$1; {s = $4+sh}; print $3, "cDNA", s, "+", 1}' 2> \
        /dev/null | sort -k1,1 -k3,3n -k4,4 | compactsga > /data/schamois/ANS_HD_UV_LD/RSEM/SGA/${filename}.mouse_cDNA.clean.Asite.sga

done
