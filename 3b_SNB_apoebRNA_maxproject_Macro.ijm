// Prompt user to select the input folder
inputDir = getDirectory("Choose input folder");

// Prompt user to select the output folder
outputDir = getDirectory("Choose output folder");

// Get a list of all image files in the input folder
fileList = getFileList(inputDir);

// Loop through all files in the folder
for (i = 0; i < fileList.length; i++) {
    filePath = inputDir + fileList[i];
    
    // Check if the file is an image (basic check for common extensions)
    if (endsWith(filePath, ".czi")) {
        
        // Open the image
        open(filePath);
        
        // Perform maximum projection
        run("Z Project...", "start=1 stop=10 projection=[Max Intensity]");
        
        // Save the projected image to the output folder with the same name
        savePath = outputDir + fileList[i];
        saveAs("Tiff", savePath);
        
        // Close all open images
        close("*");
    }
}

print("Max projection completed. Results saved to: " + outputDir);