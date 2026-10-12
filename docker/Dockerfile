FROM mambaorg/micromamba:1.5-jammy

LABEL description="pmbb-nf-toolkit-phenotyping: biobank phenotyping pipeline"
LABEL maintainer="Penn Medicine Biobank"

USER root

# Install Python dependencies via micromamba
RUN micromamba install -y -n base -c conda-forge \
    python=3.11 \
    pandas \
    numpy \
    pyyaml \
    && micromamba clean --all --yes

# Copy pipeline scripts into /opt/bin (Nextflow bin/ dir is added to PATH)
WORKDIR /opt/bin
COPY bin/ .

ENV PYTHONPATH="/opt/bin:${PYTHONPATH}"

# Make scripts executable
RUN chmod +x *.py

# Smoke tests
RUN which python \
    && python --version \
    && python -c "import pandas; print('pandas', pandas.__version__)" \
    && python -c "import numpy; print('numpy', numpy.__version__)" \
    && python -c "import yaml; print('pyyaml', yaml.__version__)"

CMD ["python", "--version"]
