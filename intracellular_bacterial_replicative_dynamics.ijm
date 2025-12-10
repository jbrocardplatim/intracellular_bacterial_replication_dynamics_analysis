//Macro Image J de mesure de prolifération bactérienne dans des cellules dont les trajectoires ont été stockées
//indépendamment via le plugin "TrackMate" dans un fichier ".xml"
//développée par Jacques Brocard, PLATIM 2020-2021, pour Florian MARRO, Evotec & CIRI
//
//Pré-requis : stack ouvert (.tif ) de plusiers canaux, les premiers = bacteria et le dernier = cellules
//et un fichier ".xml" nommé identiquement et situé dans le dossier d'origine
//
seuilb=120; //seuil detection bacteries fixe
title=getTitle();
dir=getDirectory("Image");
getDimensions(w, h, nb_channels, slices, frames);
run("Split Channels");
close();
nb_time=nSlices;
getPixelSize(unit, pixelWidth, pixelHeight);
run("Properties...", "unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1");
selectWindow("C2-" + title);
run("Smooth", "stack");
run("Subtract Background...", "rolling=25 stack");
//Ouvre le fichier ".xml" présent dans le dossier d'origine
pathfile=dir+substring(title,0,lengthOf(title)-4)+".xml";
filestring=File.openAsString(pathfile);
//Découpe chaque ligne du fichier ".xml"
rows=split(filestring, "\n");

r=1;s=0;
nLine=newArray(100); //index de ligne du début de la séquence
nSpots=newArray(100); //nb de points-temps dans la séquence
//Copie la ligne de début (nLine) et le nb de points (nSpots) pour chaque cellule trackée
while (r<rows.length-2){
	r++;
	if (substring(rows[r],lengthOf(rows[r])-12,lengthOf(rows[r])-6)=="nSpots"){
		s++;
		nLine[s]=r+1;
		nSpots[s]=substring(rows[r],lengthOf(rows[r])-4,lengthOf(rows[r])-2);
		nSpots[s]=parseInt(nSpots[s]);
	}
}
nb_cells=s;
//Pour chaque canal de bactéries, faire la mesure :
for (c=1;c<nb_channels;c++){
	rename("bacteria");
	//Fabrique le futur tableau de résultats sous la forme de nb_time+1 lignes
	tablo=newArray(nb_time+1);
	tablo[0]="Time";
	for (t=1;t<=nb_time;t++) tablo[t]=""+t;
	//Pour chaque cellule...
	for (n=1;n<=nb_cells;n++){
		//...ajouter le num de la cellule et lire les coordonnées du fichier ".xml" stockées dans "rows"
		tablo[0]=tablo[0]+"\t Cell"+n;
		verif=0;
		for (s=0;s<nSpots[n];s++){
			t=newArray(3); //tableau de coordonnées s,x,y
			t=coordon(t,rows[nLine[n]+s]); //va lire les coordonnées dans rows (fonction ci-dessous)
			posi=t[0]+1;
			//print(t[0],t[1],t[2]);
			selectWindow("bacteria");
			//Sélectionner l'image correspondant à un temps donné
			setSlice(posi);
			//Tracer un cercle de dimensions fixes
			makeOval(t[1]/pixelWidth-150,t[2]/pixelHeight-150,300,300);
			//Récupérer l'histogramme des niveaux de gris dans ce cercle...
			getHistogram(values, counts, 256);
			m=roiManager("Count");
			roiManager("Add");
			roiManager("Select",m);
			roiManager("Rename", "cell"+n+" - "+posi);
			RoiManager.setPosition(t[0]+1);
			if (t[0]>verif+1){
				for (q=verif+1;q<t[0];q++){
					tablo[q+1]=tablo[q+1]+"\t xx";	
				}
			}
			verif=t[0];
			bact=0;
			//... et compter le nb de pixels supérieur à seuilb (=bact)
			for (f=0;f<256;f++){if (values[f]>=seuilb) bact=bact+counts[f];}
			//Ajouter le résultat dans le tableau
			tablo[posi]=tablo[posi]+"\t "+bact;
		}
		if (posi<nb_time){
			for (q=t[0]+2;q<nb_time+1;q++){
				tablo[q]=tablo[q]+"\t xx";	
			}
		}
	}
	//Imprimer le tableau résultat et l'enregistrer au format texte
	for (t=0;t<=nb_time;t++) print(tablo[t]);
	pathfile=dir+substring(title,0,lengthOf(title)-4)+".txt";
	if (nb_channels>2) pathfile=dir+substring(title,0,lengthOf(title)-4)+"_c"+c+".txt";
	saveAs("Text",pathfile);
	run("Close");
	//Enregistrer puis fermer le fichier "ROI Manager"
	selectWindow("ROI Manager");
	pathfile=dir+substring(title,0,lengthOf(title)-4)+".zip";
	if (nb_channels>2) pathfile=dir+substring(title,0,lengthOf(title)-4)+"_c"+c+".zip";
	roiManager("Save",pathfile);
	run("Close");
	close("bacteria");
}

function coordon(t,r){
	i=1; lec="";
	while (lec!="t="){
		i++;
		lec=substring(r,i,i+2);
	}
	deb=i+3;
	while (lec!="x="){
		i++;
		lec=substring(r,i,i+2);
	}
	fin=i-2;
	t[0]=parseFloat(substring(r,deb,fin));
	deb=fin+5;
	while (lec!="y="){
		i++;
		lec=substring(r,i,i+2);
	}
	fin=i-2;
	t[1]=parseFloat(substring(r,deb,fin));
	deb=fin+5;
	while (lec!="z="){
		i++;
		lec=substring(r,i,i+2);
	}
	fin=i-2;
	t[2]=parseFloat(substring(r,deb,fin));
	return t;
}
