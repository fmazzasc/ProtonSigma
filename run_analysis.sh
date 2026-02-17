python3 analyse.py --output_dir "results_recal" --mom_recal
python3 analyse.py --me --output_dir "results_recal" --mom_recal
python3 ratio.py --output_dir "results_recal"


python3 analyse.py --output_dir "results_norecal"
python3 analyse.py --me --output_dir "results_norecal"
python3 ratio.py --output_dir "results_norecal"