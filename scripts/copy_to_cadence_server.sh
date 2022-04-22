#sshpass -p '30054360' scp -r $1 erfan.aghaeekiasarae@ict41702.vlsi.enel.ucalgary.ca:/nfs/atipsfs1/erfan.aghaeekiasarae/


PASSED=$1
SERVER=$2



# sshpass -p '30054360' scp $1 erfan.aghaeekiasarae@ict41702.vlsi.enel.ucalgary.ca:/nfs/atipsfs1/erfan.aghaeekiasarae/$SERVER
sshpass -p "30054360" rsync -avzhe ssh --progress $1 erfan.aghaeekiasarae@ict41704.vlsi.enel.ucalgary.ca:/nfs/atipsfs1/erfan.aghaeekiasarae/$SERVER



