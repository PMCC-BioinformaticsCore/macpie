# 1. Base image
FROM rocker/rstudio:4.3.3 AS builder
LABEL description="Docker image for macpie"

# 2. Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libz-dev \
    libomp-dev \
    cmake \
    pkg-config \
    libfontconfig1-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    libgit2-dev \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Install renv from CRAN
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# 4. Preinstall MOFA2 and Dependencies
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('MOFA2')"

# 5. Copy renv.lock and restore packages
COPY renv.lock /renv.lock
ENV renv=renv/
RUN R -e "renv::restore()"

# 6. Run MOFA2 to check installation
RUN R -e "library(MOFA2); file <- system.file(\"extdata\", \"test_data.RData\", package = \"MOFA2\"); load(file); MOFAmodel <- create_mofa(dt); MOFAmodel <- prepare_mofa(MOFAmodel); run_mofa(MOFAmodel, use_basilisk = TRUE, save_data = FALSE)"

# 7. Create final image (to reduce size)
FROM rocker/rstudio:4.3.3
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
