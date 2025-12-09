# Start from docker image from Aerts' lab
FROM aertslab/pyscenic:0.12.1

# Update and install wget
RUN apt update && apt upgrade -y && apt install -y wget

# Upgrade pip and install anndata (h5ad loading), ipykernel (jupyter), regdiffusion (to replace grnboost2)
RUN python -m pip install --upgrade pip \
    && pip install regdiffusion

# Copy the run_pyscenic.py script and make it executable
COPY run_pyscenic.py /usr/local/bin/run_pyscenic.py
RUN chmod +x /usr/local/bin/run_pyscenic.py

# Add /usr/local/bin to PATH (should already be there, but make explicit)
ENV PATH="/usr/local/bin:${PATH}"

# Create resources directory
RUN mkdir -p /opt/pyscenic_resources

# Download pySCENIC resources for human (hg38)
WORKDIR /opt/pyscenic_resources

## Feather db motif ranking files:
RUN wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
RUN wget -q "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

## MOTIF to TF annotations
RUN wget -q "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

## All TFs
RUN wget -q "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"

# Set environment variable for resource directory
ENV PYSCENIC_RESOURCES="/opt/pyscenic_resources"

# Reset working directory
WORKDIR /