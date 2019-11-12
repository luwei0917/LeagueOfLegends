# LeagueOfLegends

put all the datababse file(chromosomeA, chromosomeB, chromosomeA_to_B, chromosomeB_to_A, blastn.pkl, etc.) into datababse folder.
Everything inside folder "database, local and tracking" won't be sync.

Get result by doing:

python run_pipeline.py --AtoBTracking ~/Dropbox/genome_algorithms_comp519/project/results/GMAP_A_to_B/AtoBGMAP.tracking --BtoATracking ~/Dropbox/genome_algorithms_comp519/project/results/GMAP_B_to_A/BtoAGMAP.tracking -o test --validation dev_validation_set.tsv

python run_pipeline.py  -t C --AtoBTracking A_to_C_database/A_to_C.tracking --BtoATracking A_to_C_database/C_to_A.tracking --dbB A_to_C_database/chromosomeC --dbAtoB A_to_C_database/Minimap_chromosome_A_to_C --dbBtoA A_to_C_database/Minimap_chromosome_C_to_A  --validation A_to_C_database/test_validation_set.tsv -o test_A_to_C -n minimap_nov12 --blastnAtoB A_to_C_database/blastn_A_to_C --blastnBtoA A_to_C_database/blastn_C_to_A

Get score by doing:

python myCheck.py -e 3 ~/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_A.chromosomes/A.transcripts.fasta ~/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_B.chromosomes/B.transcripts.fasta ~/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_A.chromosomes/A.gff3 ~/Dropbox/genome_algorithms_comp519/project/Challenge_9934185_B.chromosomes/B.gff3 dev_validation_set.tsv test/GMAP_combined_nov06_post_modification_3.tsv

