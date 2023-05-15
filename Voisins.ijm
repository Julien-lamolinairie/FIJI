/*
* Macro pour compter l'aire des bulles sur plusieurs images.
*/

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".png") suffix


processFolder(input);	
// fonction pour scanner les dossiers/sous-dossiers/fichiers pour trouver tous les fichers avec le bon format
function processFolder(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder("" + input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	//Sauvegarder les résultats dans un fichier unique
	//saveAs("Results", output + "/All_Results.csv"); 

function processFile(input, output, file) {
	setBatchMode(true); // empêche l'ouverture des fenêtres d'image pendant l'exécution du script 
	// ouverture d'images via Bio-Formats
	run("Bio-Formats", "open=[" + input + "/" + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    id = getImageID(); // obtenir l'id de l'image de base
	run("Duplicate...", " "); // dupliquer cette image pour travailler sur sa copie
	original=getTitle();
    
	getDimensions(width, height, channels, slices, frames);
	run("Analyze Particles...", " show=Masks display exclude clear");
	run("Invert LUT");
	rename(original+"-1");
	original=getTitle();
	run("Duplicate...", "title=[NbHood_"+original+"]");
	neighborhood=getTitle();

	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear record");
	//define variables
	initialParticles=nResults;
	X=newArray(nResults);
	Y=newArray(nResults);
	neighborArray=newArray(nResults);
	neighbors=0;
	mostNeighbors=0;
	//retveive particle coordinates
	for(l=0; l<initialParticles; l++) {
		X[l]=getResult("XStart", l);
		Y[l]=getResult("YStart", l);
	}
	//prepare selector image
	setForegroundColor(255, 255, 255);
	setBackgroundColor(0, 0, 0);
	run("Set Measurements...", " centroid redirect=None decimal=3");
	run("Wand Tool...", "mode=8-connected tolerance=0");
	run("Options...", "iterations=1 count=1 black edm=Overwrite do=Nothing");
		

	for(hood=0; hood<initialParticles; hood++) {
			//create selector neighborhood
			selectWindow(neighborhood);
			run("Select None");
			run("Duplicate...", "title=[Selector_"+original+"]");
			selector=getTitle();
			doWand(X[hood], Y[hood], 0, "8-connected");
			//print(hood + "(" + X[hood]+"/"+ Y[hood] + ")");
			run("Enlarge...", "enlarge="+20 + " pixel");
			run("Fill");
			run("Select None");
			doWand(X[hood], Y[hood], 0, "8-connected");

			selectWindow(neighborhood);

			run("Restore Selection");
			run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing clear record");
			neighbors = nResults-1;
			neighborArray[hood]=neighbors;
			if(neighbors>mostNeighbors) {
				mostNeighbors=neighbors;
			}
			close(selector);
		}
		
	selectWindow(original);
	run("Duplicate...", "title=[P-NbHood_"+20+"_"+original+"]");
	particles=getTitle();
	selectWindow(particles);

	for(mark=0; mark<initialParticles; mark++) {
		markValue=neighborArray[mark];
		if(markValue==0) {
			doWand(X[mark],Y[mark], 0, "8-connected");
			Roi.setStrokeColor(0);
			run("Add Selection...");
			run("Select None");
		}
		setForegroundColor(markValue, markValue, markValue);
		floodFill(X[mark],Y[mark], "8-connected");
		
	}

	
	run("Select None");		
	run("glasbey");
	setBatchMode("show");
		
	//visually eliminate edge particles (but count them as neighbors)

	selectWindow(original);
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing exclude clear record");
	visibleParticleNumber = nResults;
	
	

	//create distribution plot

	selectWindow(original);
	run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 pixel show=Masks exclude clear record add");
	roiManager("Show None");
	masksWithoutEdgeParticles = getTitle();
	run("Select All");
	run("Copy");
	selectWindow(particles);
	setPasteMode("Transparent-white");
	run("Paste");
	setPasteMode("Copy");
	colorIndex = 0;
	neighborList = newArray(mostNeighbors+1);
		for(roi=0; roi<roiManager("Count"); roi++) {
			selectWindow(particles);
			roiManager("Select", roi);
			getStatistics(areaNotNeeded, colorIndex);
			neighborList[colorIndex] = neighborList[colorIndex] + 1;
		}
		particleCount = roiManager("Count");
	} else {
		neighborList = newArray(mostNeighbors+1);
		Array.fill(neighborList, 0);
		for(num=0; num<initialParticles; num++) {
			nextNeighbor = neighborArray[num];
			if(nextNeighbor>0) {
				neighborList[nextNeighbor] += 1;
			} else {
					
neighborList[0] += 1;

			}
		}
		particleCount = initialParticles;
	}
		
		
	Plot.create("Distribution: " + particles, "neighbors", "count", neighborList);
	Plot.addText("particles (total) = " + particleCount, 0.01, 0.1);
	setBatchMode("show");
}
	

	
//Calibration Bar
if(calibrationbar==true) {
	stepsize=floor(256/mostNeighbors);
	newImage("Calibration_"+original, "8-bit Black", (stepsize*mostNeighbors+stepsize), 30, 1);
	w=getWidth();
	step=0;
	for(c=0; c<=mostNeighbors+1; c++) {
		makeRectangle(step, 0, step+stepsize, 30);
		setForegroundColor(c, c, c);
		run("Fill");
		step=step+stepsize;
	}
	run("Select None");
	run("glasbey");
	run("RGB Color");
	setForegroundColor(255, 255, 255);
	setBackgroundColor(0, 0, 0);
	run("Canvas Size...", "width="+w+" height=50 position=Top-Center");
	if(mostNeighbors>9) { 
		offset=15;
	} else {
		offset=10;
	}
	drawString("0", 2, 48);
	drawString(mostNeighbors, w-offset, 48);
}
setBatchMode("show");


}




	save(output + "/Binary_OUTPUT_" + file);
    run("Analyze Particles...", "size=10-infinity exclude add");
    selectImage(id); // activer l'image de base
    roiManager("Show All with labels"); // overlay ROIs
	roiManager("Deselect");
	roiManager("Measure"); // measure on original image
	saveAs("Results", output + "/"+file+"_All_Results.csv"); 
	run("Clear Results");
	// save ROIs for current image
	roiManager("Deselect");
	roiManager("Save", output+ "/" + file + "_ROI.zip"); // saves Rois zip file
	roiManager("Deselect");
	roiManager("Delete"); // clear ROI Manager for next image
}


