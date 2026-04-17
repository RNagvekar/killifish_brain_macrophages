# Define the function to reformat the text
reformat_text_file <- function(input_file, output_file) {
  # Read the content of the text file
  text_content <- readLines(input_file, warn = FALSE)
  
  # Collapse the lines into a single string and remove spaces and line breaks
  cleaned_text <- gsub("[[:space:]]+", "", paste(text_content, collapse = ""))
  
  # Split the cleaned text into chunks of 60 characters
  formatted_lines <- strsplit(cleaned_text, "(?<=.{60})", perl = TRUE)[[1]]
  
  # Write the formatted lines to the output file
  writeLines(formatted_lines, output_file)
}

# Specify the input and output file paths
input_file_path <- "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/Figures/Raw+processed_data/Genotyping/oScarletDNA/20250901_oScarlet_insertion_translated_oScarletinframe.txt"
output_file_path <- "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/Figures/Raw+processed_data/Genotyping/oScarletDNA/20250901_oScarlet_insertion_translated_oScarletinframe_60perline.txt"

# Call the function to reformat the text file
reformat_text_file(input_file_path, output_file_path)

# Print a message indicating completion
cat("Text reformatting complete. Output saved to:", output_file_path, "/n")