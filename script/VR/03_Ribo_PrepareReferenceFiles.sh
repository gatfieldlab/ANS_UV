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

# note about the transcriptCodonsNorm.sga, found in a RUST analysis script:
# # 1. Generate target feature file. Here it is the codons, but this can
# #    be changed to any feature that can be mapped to a single base.
/path_to/makeCodonFile.py -g /data/databases/mouse/gtf/Mmusculus.GRCm39.111.gtf \
                -f /data/databases/mouse/fasta/Mmusculus.GRCm39.111.cdna.ensembl.fa \
                > ./config/Mmusculus.GRCm39.111.transcriptCodons.sga
/path_to/calculateCodonFrequency.py -n -i config/Mmusculus.GRCm39.111.transcriptCodons.sga > \
                          ./config/Mmusculus.GRCm39.111.transcript.normFactor
 awk 'BEGIN{FS=OFS="\t"; while((getline < "config/Mmusculus.GRCm39.111.transcript.normFactor")>0){nf[$1][$2]=$3}}{print $1,$2,$3,$4,nf[$1][$2]}' config/Mmusculus.GRCm39.111.transcriptCodons.sga>   config/Mmusculus.GRCm39.111.transcriptCodonsNorm.sga
