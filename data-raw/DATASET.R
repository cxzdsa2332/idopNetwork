# Load raw data from .csv file
gut_microbe <- read.csv("data-raw/gut_microbe.csv")
mustard_microbe <- read.csv("data-raw/mustard_microbe.csv", row.names = 1, check.names = F)
# Apply preprocessing...
# Save the cleaned data in the required R package location
usethis::use_data(gut_microbe, overwrite = TRUE)
usethis::use_data(mustard_microbe, overwrite = TRUE)
