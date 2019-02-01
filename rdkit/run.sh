apt-get install wget
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh
source ~/.bashrc
conda install -c rdkit rdkit 
cd /home/ && git clone https://github.com/n-yoshikawa/openbabel-paper.git
cd openbabel-paper/rdkit/
{ time python convert-rdkit.py; } 2> time-log
