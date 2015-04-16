storepath=/hdfs/store/user/khurana/
dirname=$1
filename=dirname.txt
rm $filename
fullpath=$storepath/$dirname
ls -1 $fullpath >& $filename
for which in `less $filename`
do 
    files=`find  $fullpath/$which  -name "*.root" | gawk '{ORS=" "}{print $1}'`
    hadd Merged_${which}.root $files
done 
