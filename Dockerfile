# Force x86_64 so INLA pre-built binaries are available (stable repo has no aarch64 builds)
FROM --platform=linux/amd64 rocker/r-ver:4.5.1

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libabsl-dev \
    cmake \
    wget \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Use Posit Package Manager for pre-compiled binaries (Ubuntu Noble / R 4.5)
ENV CRAN="https://packagemanager.posit.co/cran/__linux__/noble/latest"

# Install R packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'magrittr', \
    'splines', \
    'ggplot2', \
    'BiocManager', \
    'fixest', \
    'glmmTMB', \
    'MuMIn' \
  ), repos = Sys.getenv('CRAN'), Ncpus = parallel::detectCores())" \
  && R -e "stopifnot(requireNamespace('tidyverse'), requireNamespace('fixest'), requireNamespace('glmmTMB'))"

# Install fmesher explicitly (from CRAN) BEFORE INLA
RUN R -e "install.packages('fmesher', repos=Sys.getenv('CRAN'), Ncpus=parallel::detectCores())" \
  && R -e "stopifnot(requireNamespace('fmesher'))"

# Install INLA from its own repository
# INLA is not on CRAN so must be installed this way
# CRAN repo is included so dependency 'fmesher' can be resolved
RUN R -e "options(timeout=600); install.packages( \
    'INLA', \
    repos = c(INLA = 'https://inla.r-inla-download.org/R/testing', CRAN = Sys.getenv('CRAN')), \
    dependencies = TRUE)" \
  && R -e "stopifnot(requireNamespace('INLA'))"

# Fix INLA binary permissions so Singularity can execute them
RUN R -e "path <- system.file('bin', package='INLA'); if (nchar(path)) system(paste('chmod -R 755', path)) else message('INLA bin/ not found — skipping')"

# Verify INLA installed correctly
RUN R -e "library(INLA); inla.version()"

RUN apt-get update && apt-get install -y \
    python3 \
    python3-venv \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Python virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Python packages
COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt jupyterlab

# Run environment tests
COPY tests/ /tmp/tests/
RUN Rscript /tmp/tests/test_r_packages.R
RUN python3 /tmp/tests/test_python_packages.py

WORKDIR /project
RUN mkdir -p /project/data /project/scripts /project/outputs

EXPOSE 8888
CMD ["R"]
