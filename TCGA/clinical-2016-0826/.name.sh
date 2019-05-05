awk '{
        grep="str=`grep \"admin:disease_code xsd_ver=\" "$1"/"$2"`";
	str="str1=${str#*>}";
	str1="str2=${str1%<*}";
	echo="echo $str2";
	mkdir="mkdir $str2";
	file=$1"/"$2;
	mv=sprintf("mv -f %s $str2/",file);
	rm=sprintf("rm -rf %s",$1);
	cmd=sprintf("%s\n%s\n%s\n%s\n%s\n%s",grep,str,str1,echo,mkdir,mv);
	print cmd;
	system(cmd);
}' gdc0826.txt
