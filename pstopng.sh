folder_ps=Analyse/PS1
folder_png=png
mkdir -p $folder_png ;
for i in $(ls $folder_ps/*ps) ; do  name_png=$(echo $i | cut -f 1 -d '.').png ; echo "$i to $name_png" ; convert -density 10 $i $name_png  ; mv $name_png $folder_png ; done
echo "test"
