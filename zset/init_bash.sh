path=$(cd "$(dirname $0)";pwd)
p_path=$(cd $path ;cd ../;pwd)
echo "# >>> zhangxin alias:$(date "+%Y-%m-%d-%H-%M-%S") >>>" >> ~/.bashrc
for i in $(cd $p_path;ls zbash;)
do
	echo "alias ${i}='bash ${p_path}/zbash/${i}'" >> ~/.bashrc
done
echo "# <<< zhangxin alias:$(date "+%Y-%m-%d-%H-%M-%S") <<<" >> ~/.bashrc
