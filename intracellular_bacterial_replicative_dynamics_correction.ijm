// Compter les pixels >= 140 dans ROI defini

	{

    seuilb = 120;  //seuil detection bacteries fixe

    //Récupérer l'histogramme des niveaux de gris dans le ROI defini
    getHistogram(values, counts, 256);
	//Et compter le nb de pixels supérieur à seuilb (=bact)
    bact = 0;
    for (f = 0; f < values.length; f++) {
        if (values[f] >= seuilb)
            bact = bact + counts[f];
    }
    
    // Résultats
    print("Nombre de pixels >= " + seuilb + " dans le ROI : " + bact);
