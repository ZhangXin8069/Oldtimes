home='/content/drive/MyDrive/'
cd $home
pip install --upgrade pip
sudo rm -rf /content/drive/MyDrive/zhangxin
git clone https://github.com/zx2810/zhangxin.git
cp /content/drive/MyDrive/zhangxin/* /content/drive/MyDrive/
sudo rm -rf /content/drive/MyDrive/zhangxin
pip install -r zhangxin_requirements.txt
sudo apt update
sudp apt install git cmake g++ mpich gfortran libhdf5-dev libpng-dev libeigen3-dev libmpich-dev libnetcdf-dev libtbb-dev libgles2-mesa-dev
pip install -r openmc_requirements.txt
git clone https://github.com/openmc-dev/openmc.git
cd openmc
sudo make install
pip install . && cmake .

