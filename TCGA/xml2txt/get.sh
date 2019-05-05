fi=`ls TCGA/TCGA/clinical-2016-0826`
for e in ${fi[@]}
do
	cd TCGA/TCGA/clinical-2016-0826/$e
	file=`ls *.xml`
	for ee in ${file[@]}
	do
		python TCGA/TCGA/xml2txt/xml2txt.py -i $ee -o TCGA/TCGA/xml2txt/$e/$ee.txt
		if [ "$?" = 0 ]
		then
			echo "$e/$ee successfully"
		else
			echo "$e/$ee failed " >>TCGA/TCGA/xml2txt/error.log
		fi
	done
	echo "$e complete"
done
