/*
 * Nataliya Trushina, 2020
 * Fiji stress granule image analysis: macro processing multiple images in a folder
 * Requirements: nd2 (or tiff) files with image series containing one image before acquisition and a series of images acquired every second starting one second after photoactivation
 * Output: csv files with measured intensities and archives with ROIs assigned during analysis 
 */

// Check that there are not any files with the chosen extension in subfolders or make sure that the analysis is also required for images in subfolders

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".nd2") suffix

run("Set Measurements...", "integrated redirect=None decimal=3");
roiManager("reset");
processFolder(input);

// function to scan folders/subfolders/files to find files with selected suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	print("Processing: " + input + File.separator + file);
	open(input + File.separator + file);
	selectWindow(file +" - C=0"); // choose the channel for acquisition, might need to be changed when using PA-mCherry for example
	Stack.setFrame(2); // select second frame (first frame acquired after photoactivation)
	run("Enhance Contrast", "saturated=0.35"); // adjust brigthness for second frame automatically 
	name = File.nameWithoutExtension();
	id=getImageID();
	// the instructions appear in the new windows
	waitForUser("Draw ROI and check, close file if ROI is inappropriate");
	if (isOpen(id) == false) {
	 	f = File.open(output + File.separator + name + "_inappr_ROI" + ".txt");
	 	File.close(f);
	}
	else {
	// all the analysis steps if ROI is selected
	roiManager("Add");
	roiManager("Multi Measure");
	saveAs("Results", output + File.separator + name + "_RawIntDen1.csv");
	roiManager("Save", output + File.separator + name + ".zip");
	roiManager("Delete");
	close();
	}
}
print("The analysis is complete.")