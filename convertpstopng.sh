#ffmpeg -start_number 126 -i img%03d.png -pix_fmt yuv420p out.mp4

density=20
for INP in *.ps
do
	newname=$(echo "$INP" | cut -d'.' -f1)
	newname=$(echo "$newname.png")
	convert -density $density -geometry 100% $INP $newname
	echo "convert $INP to $newname complete"
done
