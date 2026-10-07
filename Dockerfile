# NeuroGut-MetaSeq Reproducible Container Environment
FROM python:3.11-slim

LABEL maintainer="NeuroGut-MetaSeq Team"
LABEL description="Containerized environment for reproducible RNA-seq meta-analysis of microglia"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install essential build and runtime libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgomp1 \
    make \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Install Python dependencies
COPY requirements.txt /workspace/
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy source repository
COPY . /workspace/

# Default entrypoint: verify and build paper
CMD ["make", "paper"]
