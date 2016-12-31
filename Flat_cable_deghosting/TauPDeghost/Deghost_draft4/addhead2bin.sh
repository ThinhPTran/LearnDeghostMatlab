if [ "$#" -ne 5 ]; then
  echo "Usage: adhead2bin binfile headerfile Nt SuFile SegyFile"
  exit 1
fi

echo "bin file: " $1
echo "header file: " $2
echo "Nt: " $3
echo "Su file: " $4
echo "Segy file: " $5
 
grep -v [a-zA-Z] $2 > tmp.ascii
a2b <tmp.ascii n1=2 > tmp.ascii.bin
suaddhead ns=$3 < $1 | sushw key=dt a=4000 > tmp.su
sushw <tmp.su key=offset,cdp infile=tmp.ascii.bin > $4
rm -f tmp.ascii.bin tmp.su tmp.ascii 
surange <$4
segyhdrs <$4 hfile=header bfile=binary
segywrite <$4 hfile=header bfile=binary tape=$5
rm -f header binary
