if [ -f /unsup/R/bin/R ]
then
	mkdir output

    /unsup/R/bin/Rscript inference.r $*

	tar -zcvf output/output-$1.tgz2 output/*.txt
	rm output/*.txt

	exit 0
else 
     exit 1
fi
