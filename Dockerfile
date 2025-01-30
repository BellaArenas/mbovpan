# Use an official Miniconda image as the base
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /mbovpan

# Create a Conda environment with Python 3.10 and Mamba installed
RUN conda create -n mbovpan python=3.10 mamba -c conda-forge -y

# Install required packages in the mbovpan environment using mamba for faster installation
RUN conda run -n mbovpan mamba install -c conda-forge -c bioconda nextflow=22.10.6 pandas panaroo sra-tools -y

# Clean up to reduce image size
RUN conda clean -afy

# Clone the mbovpan repository (if needed)
# Note: Adjust the repository URL as necessary.
RUN git clone https://github.com/bellaarenas/mbovpan.git /mbovpan

# Set the PATH to include the mbovpan environment's bin directory
ENV PATH /opt/conda/envs/mbovpan/bin:$PATH

# Create mbovis_input directory and download sequences
RUN mkdir /mbovpan/mbovis_input && \
    cd /mbovpan/mbovis_input && \
    fasterq-dump --verbose --split-3 SRR10482974 && \
    fasterq-dump --verbose --split-3 SRR10482944 && \
    fasterq-dump --verbose --split-3 ERR11893527 && \
    fasterq-dump --verbose --split-3 SRR23174187 && \
    fasterq-dump --verbose --split-3 SRR23174144 && \
    cd /mbovpan

# Set the default command to run when starting the container
CMD ["bash"]