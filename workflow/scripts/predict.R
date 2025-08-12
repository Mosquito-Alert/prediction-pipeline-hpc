# Load libraries
library(brms)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3){
    stop("Please provide exactly three arguments: model_rds_file input_features.csv output.csv")
}
model_rds_file <- args[1]
input_csv <- args[2]
output_csv <- args[3]

# Load the saved model
model <- readRDS(model_rds_file)

# Load data for prediction
data <- read.csv(input_csv)

# Create a new column 'SE' and set to 1
data$SE <- 1

# Predict on new data
predictions <- posterior_predict(
    model,
    newdata = data,
    allow_new_levels = TRUE,
    re_formula =  NA,
    ndraws = 1000
)

pred_mean <- colMeans(predictions)

# Create output dataframe
output_df <- data.frame(
    id = data$id,
    prediction = pred_mean
)

# Write to CSV
write.csv(output_df, output_csv, row.names = FALSE)