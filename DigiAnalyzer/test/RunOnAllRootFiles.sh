i=0
for which in `less filelist.txt`; do
    rm digianalyzer_${i}_cfg.py
    sed 's/INPUTFILENAME/'${which}'/' <digianalyzer_cfg.py > digianalyzer_${i}_cfg.py
    cmsRun digianalyzer_${i}_cfg.py
    i=$(($i+1))
done