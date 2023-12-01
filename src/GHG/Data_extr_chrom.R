# Load required libraries
# install.packages("pdftools")
# install.packages("stringr")
# install.packages("dplyr")
library(pdftools)
library(stringr)
library(dplyr)
# devtools::install_github("ropensci/tabulizer")
# install.packages("rJava")
# Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-18/")
library(rJava)
# install.packages("tabulizer")
library(tabulizer)
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

# Function to extract data from each PDF file and create a data frame
extract_data_from_pdf <- function(pdf_file_path) {
  pdf_name <- str_extract(basename(pdf_file_path), "\\d{8}")
  
  # Extract data from the Front Signal Results table
  front_table <- extract_tables(pdf_file_path, pages = 1)
  if (length(front_table) > 0) {
    front_table <- front_table[[1]]
    if (nrow(front_table) >= 2) {
      col_names <- front_table[1, ]
      col_names <- str_to_lower(col_names)
      col_names <- str_replace_all(col_names, "[[:punct:]]", "")
      col_names <- str_trim(col_names)
      
      front_table <- front_table[-1, ]
      col_index <- which(col_names %in% c("name", "retention time", "area", "concentration"))
      
      if (length(col_index) == 4) {
        front_table <- front_table[, col_index]
        colnames(front_table) <- c("Name", "Retention Time", "Area", "Concentration")
        
        CH4_ppm <- as.numeric(front_table[front_table$Name == "CH4", "Concentration"])
        CO2_ppm <- as.numeric(front_table[front_table$Name == "CO2", "Concentration"])
        N2O_ppm <- as.numeric(front_table[front_table$Name == "N2O", "Concentration"])
      } else {
        CH4_ppm <- NA
        CO2_ppm <- NA
        N2O_ppm <- NA
      }
    } else {
      CH4_ppm <- NA
      CO2_ppm <- NA
      N2O_ppm <- NA
    }
  } else {
    CH4_ppm <- NA
    CO2_ppm <- NA
    N2O_ppm <- NA
  }
  
  # Combine the data into a data frame
  data <- data.frame(PDF_Name = pdf_name,
                     CH4_ppm = CH4_ppm,
                     CO2_ppm = CO2_ppm,
                     N2O_ppm = N2O_ppm)
  return(data)
}

# Set the folder path where your PDF files are located.
pdf_folder <- "C:/Users/SECHEVERRIA/IRTA/Gas results/Chromatography/Soporte_JN/Chromat_reports"

# List all the PDF files in the folder and apply the function to each file.
pdf_files <- list.files(pdf_folder, pattern = "*.pdf", full.names = TRUE)
all_data <- lapply(pdf_files, extract_data_from_pdf)

# Combine the individual data frames into a single data frame.
combined_data <- do.call(rbind, all_data)

# Print the result
print(combined_data)