content="~/content/"
zbash="${content}zbash/"
mkdir -p zbash
mv 'z*\.sh' $zbash
cp ~/.bashrc $zbash	
echo "hao" >> ~/.bashrc
source ~/.bashrc
