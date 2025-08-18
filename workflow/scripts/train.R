# Use snakemake S4 object
features_file <- snakemake@input[[1]]
output_file   <- snakemake@output[[1]]

# Load necessary libraries
library(brms)
library(cmdstanr)
library(readr)
library(dplyr)
library(lubridate)

# Set MCMC parameters
iteret <- 5000
wup <- 2000
nchains <- snakemake@threads
threads_per_chain <- 1

# Load data
features <- read_csv(features_file) %>%
    mutate(
        year = year(as.Date(date)),
        has_bite_reports = n_bite_reports > 0
    )

# features_clean <- features[!is.na(features$min_t2m_21d), ]
features_clean <- features[!is.na(features$SE), ]

model <- brm(
    has_bite_reports ~ poly(min_t2m_21d, 2) + 
        perc_agricultural + perc_other  + perc_discont_urban_fabric +
        (1 | h3_index) + (1 | year) +
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