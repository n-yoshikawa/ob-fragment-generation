python calc_rmsd.py ../OB-new/platinum-ob-new.sdf > rmsd-ob-new.csv
python calc_rmsd.py ../OB-org/platinum-ob-org.sdf > rmsd-ob-org.csv
python calc_rmsd.py ../rdkit/platinum-rdkit-etkdg.sdf > rmsd_rdkit.csv
./eval ../data/platinum_dataset_2017_01.sdf ../rdkit/platinum-rdkit-etkdg.sdf > bond-angle-rdkit.csv 2>/dev/null
./eval ../data/platinum_dataset_2017_01.sdf ../OB-new/platinum-ob-new.sdf > bond-angle-ob-new.csv 2>/dev/null
./eval ../data/platinum_dataset_2017_01.sdf ../OB-org/platinum-ob-org.sdf > bond-angle-ob-org.csv 2>/dev/null
python plot-result.py rmsd-ob-new.csv bond-angle-ob-new.csv plot-rmsd-ob-new.pdf "New Open Babel"
python plot-result.py rmsd-ob-org.csv bond-angle-ob-org.csv plot-rmsd-ob-org.pdf "Original Open Babel"
python plot-result.py rmsd_rdkit.csv bond-angle-rdkit.csv plot-rmsd-rdkit.pdf "RDKit"
