<<<<<<< HEAD
#!bin/sh
export LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | tr ":" "\n" | grep -v "stubs" | tr "\n" ":"`
python ./esmfold/esmfold.py -i ../Example/7stb_C.fasta -o ../Example/structure_data --chunk-size 128
# python extract.py -i ./DIPS.fasta -o ./ldq_fasta/
#--chunk-size 128
#--cpu-offload
=======
#!bin/sh
export LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | tr ":" "\n" | grep -v "stubs" | tr "\n" ":"`
python ./esmfold/esmfold.py -i ../Example/7stb_C.fasta -o ../Example/structure_data --chunk-size 128
# python extract.py -i ./DIPS.fasta -o ./ldq_fasta/
#--chunk-size 128
#--cpu-offload
>>>>>>> b517712a90b961091f509de7b63373e7c68c0d25
