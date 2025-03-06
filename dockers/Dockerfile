FROM rocker/rstudio:4.3.3
LABEL description="Docker image for macpie"

ENV R_MAKEVARS="-std=c++17"

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


# 2. Install renv from CRAN
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# 2.1. Force R to use C++17 in Makeconf
#    - We override all relevant variables so we don't end up with -std=gnu++11
RUN echo "CXX11 = g++"                         >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXX11FLAGS = -std=c++17 -fPIC"        >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXX14 = g++"                         >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXX14FLAGS = -std=c++17 -fPIC"        >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXX17 = g++"                         >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXX17FLAGS = -std=c++17 -fPIC"        >> /usr/local/lib/R/etc/Makeconf && \
    echo "CXXFLAGS = -std=c++17 -fPIC"          >> /usr/local/lib/R/etc/Makeconf


# Install MOFA2 and Dependencies
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('MOFA2')"

# 3. Copy your renv.lock file into the container
COPY renv.lock /renv.lock
ENV renv=renv/

# 4. Restore R packages from renv.lock
RUN R -e "renv::restore()"

# Run MOFA2 to check installation
RUN R -e "library(MOFA2); file <- system.file(\"extdata\", \"test_data.RData\", package = \"MOFA2\"); load(file); MOFAmodel <- create_mofa(dt); MOFAmodel <- prepare_mofa(MOFAmodel); run_mofa(MOFAmodel, use_basilisk = TRUE, save_data = FALSE)"

