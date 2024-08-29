FROM continuumio/miniconda3

# Set working directory
WORKDIR /src

# Install RDKit from conda-forge channel
RUN conda install conda-forge::rdkit

# Copy requirements file
COPY requirements.txt requirements.txt

# Upgrade pip and install Python dependencies
RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Expose port
EXPOSE 8000

# Entry point for container
ENTRYPOINT [ "uvicorn" ]

# Run app
CMD [ "main:app", "--host", "0.0.0.0", "--port", "8000" ]