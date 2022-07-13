normalized_counts_file=$1
annotationfile=$2
normalized_counts_annotation_file=$3
head -n 1 $normalized_counts_file > $normalized_counts_annotation_file
sort -u -k1,1 $annotationfile >> tmp.txt
mv tmp.txt $annotationfile
sed -i "s/\tNA\tNA\t/\t/" $annotationfile
cat $annotationfile | while read var1 var2 var3 var4 var5;do grep $var1 $normalized_counts_file | sed "s/$var1/$var2/" >> $normalized_counts_annotation_file ;done
