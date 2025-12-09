# Start from docker image from Aerts' lab
FROM aertslab/pyscenic:0.12.1

# Update
RUN apt update && apt upgrade -y

# Upgrade pip and install anndata (h5ad loading), ipykernel (jupyter), regdiffusion (to replace grnboost2)
RUN python -m pip install --upgrade pip \
    && pip install regdiffusion