density=20
for INP in *.ps
do
	newname=$(echo "$INP" | cut -d'.' -f1)
	newname=$(echo "$newname.png")
	convert -density $density -geometry 100% $INP $newname
	echo "convert $INP to $newname complete"
done
