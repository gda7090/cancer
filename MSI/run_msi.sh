
echo Starting Analysis of $PFX > $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;

echo "Making mpileups" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM | awk '{if($4 != 0) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 

echo "Varscan Readcounts start" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
wait

echo "MSI Analyzer start" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;

msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt

echo "MSI calls start" >> $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;

msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt
rm -f $SAVEPATH/$PFX/$PFX.mpileup
echo Completed Analysis of $PFX>> $SAVEPATH/$PFX/"$PFX"_msi_log;
date +"%D %H:%M" >> $SAVEPATH/$PFX/"$PFX"_msi_log;

exit

#msi count_msi_samples $MSI_BASELINE $SAVEPATH -o $SAVEPATH/Combined_MSI.txt



