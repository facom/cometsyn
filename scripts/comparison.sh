files=$(ls -m output/version2/*NOV_13*.dat)
python comet-properties.py BASEDIR=runs/version2/run-nov13-rock ExpansionComparison "files='$files'"
