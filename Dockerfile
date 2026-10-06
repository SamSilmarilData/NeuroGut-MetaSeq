# NeuroGut-MetaSeq Reproducible Environment
FROM rocker/r-ver:4.3.2

LABEL maintainer="NeuroGut-MetaSeq Team"
LABEL description="Containerized environment for reproducible RNA-seq meta-analysis of microglia"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    make \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Install R packages
COPY install_packages.R /workspace/
RUN Rscript /workspace/install_packages.R

# Install Python packages
COPY requirements.txt /workspace/
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r requirements.txt

# Copy source code and config
COPY . /workspace/

CMD ["make", "paper"]
