source venv/bin/activate

#ENV Variables
DATA="/home/isaac/dev/sfu/cmpt415/hest_data"
PROJECT=$(pwd)
echo "Setting up Cedar Filesystem"
sshfs -v cedar:/home/imv/scratch data
