path=$(cd "$(dirname $0)";pwd)
p_path=$(cd $path ;cd ../;pwd)
log_path=${path}/alias_log.txt
rm $log_path
for i in $(cd $p_path;ls zbash;)
do
	echo "alias ${i}='bash ${p_path}/zbash/${i}'" >> $log_path
done
bash $log_pathe
