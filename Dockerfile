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
	&& pip install --no-cache-dir uv \
	&& rm -rf /var/lib/apt/lists/*

# Resources - download all in single RUN command (before switching to non-root user)
RUN mkdir -p /opt/pyscenic_resources && \
	cd /opt/pyscenic_resources && \
	wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" && \
	wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" && \
	wget -q "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" && \
	wget -q "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"

# Copy project metadata and install deps into an isolated venv (non-root friendly)
WORKDIR /app
COPY pyproject.toml README.md run_pyscenic.py ./
RUN python -m venv /opt/venv \
	&& /opt/venv/bin/pip install --no-cache-dir uv \
	&& /opt/venv/bin/uv pip install --no-cache .

# Use non-root user for HPC/Singularity compatibility
RUN useradd -m -u 1000 -s /bin/bash pyscenic \
	&& chown -R pyscenic:pyscenic /app /opt/venv /opt/pyscenic_resources

USER pyscenic
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
ENV PYSCENIC_RESOURCES="/opt/pyscenic_resources"

WORKDIR /app
ENTRYPOINT ["run-pyscenic"]