macro_code <- '
inputDirectory = "/home/bacteraia/transwell kd hif/LZ/SC/1/";
outputDirectory = "/home/bacteraia/transwell kd hif/LZ/SC/1/";

// Get list of files
list = getFileList(inputDirectory);

for (i = 0; i < list.length; i++) {
  open(inputDirectory + list[i]);
  
  // Convert image to grayscale if needed
  run("8-bit");
  
  // Adjust threshold (try different values for best results)
  setAutoThreshold("Default");
  
  // Invert LUT if necessary
  // run("Invert");
  
  // Perform Watershed (optional, useful for splitting merged cells)
  run("Watershed");
  
  // Analyze particles (set size and circularity as needed)
  run("Analyze Particles...", "size=300-Infinity circularity=0.30-1.00 show=Overlay display exclude summarize");
  
  // Save the results
  saveAs("Results", outputDirectory + list[i] + "_results.csv");
  
  // Close image
  close();
}
'

macro_path <- "/home/bacteraia/transwell kd hif/counting.ijm"
writelines(macro_code, con = macro_path)
