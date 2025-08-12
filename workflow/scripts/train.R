# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
    stop("Please provide the path to the features CSV file as the first argument.")
}
features_file <- args[1]
output_file <- args[2]

# Load necessary libraries
library(brms)
library(cmdstanr)
library(readr)
library(dplyr)
library(lubridate)

# Set MCMC parameters
iteret <- 5000
wup <- 2000
nchains <- 4
threads_per_chain <- 1

# Load data
features <- read_csv(features_file) %>%
    mutate(year = year(as.Date(date)))

features_clean <- features[!is.na(features$min_t2m_21d), ]
features_clean <- features[!is.na(features$SE), ]

model <- brm(
    has_bite_reports ~ poly(min_t2m_21d, 2) + 
        perc_agricultural + perc_other  + perc_discont_urban_fabric +
        (1 | id) + (1 | year) + 
        offset(log(SE)),
    data = features_clean,
    prior = set_prior("cauchy(0,2.5)", class="b"),
    family = bernoulli(link = "logit"),
    iter = iteret,
    warmup = wup,
    chains = nchains,
    cores = nchains,
    backend = "cmdstanr",
    threads = threading(threads_per_chain),
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99))

saveRDS(model, file = output_file)