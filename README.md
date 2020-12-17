# pairwise2_aligner_GUI

venv
    python 3.8.6
    pip install numpy==1.19.1
    pip install pandas==1.1.1
    pip install biopython==1.77
    pip install pyinstaller==4.0
    pip install xlrd==1.2.0
    pip install openpyxl==3.0.5
     
how to make .exe file

    pyinstaller --onefile --noconsol --icon=dna.ico aligner.py
    pyinstaller --onefile --noconsol --icon=dna.ico --add-data example.PNG;. aligner.py

    then check dist folder

[excute file link](./aligner.exe)
