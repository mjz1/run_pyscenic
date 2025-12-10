# Base image: minimal Python 3.10
FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive \
	PYTHONDONTWRITEBYTECODE=1 \
	PYTHONUNBUFFERED=1

# System deps for scientific Python and building wheels, plus uv installation
RUN apt-get update && apt-get install -y --no-install-recommends \
	build-essential \
	curl \
	wget \
	git \
	ca-certificates \
	libhdf5-dev \
	libxml2-dev \
	libxslt1-dev \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libssl-dev \
	libffi-dev \
	liblapack-dev \
	libblas-dev \
	gfortran \
	pkg-config \
	&& rm -rf /var/lib/apt/lists/* \
	&& curl -LsSf https://astral.sh/uv/install.sh | sh

# Resources - download all in single RUN command (before switching to non-root user)
RUN mkdir -p /pyscenic/resources && \
	cd /pyscenic/resources && \
	wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" && \
	wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" && \
	wget -q "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" && \
	wget -q "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"

# Copy the project into the image
COPY . /pyscenic

# Sync the project into a new environment, asserting the lockfile is up to date
WORKDIR /pyscenic
RUN uv sync

# Use non-root user for HPC/Singularity compatibility
RUN useradd -m -u 1000 -s /bin/bash pyscenic \
	&& chown -R pyscenic:pyscenic /pyscenic
USER pyscenic

ENV VIRTUAL_ENV=/pyscenic/.venv
ENV PATH="/pyscenic/.venv/bin:${PATH}"
ENV PYSCENIC_RESOURCES="/pyscenic/resources"

ENTRYPOINT ["run-pyscenic"]