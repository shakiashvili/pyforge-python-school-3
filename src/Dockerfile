# Use Minicoda image  
FROM continuumio/miniconda3
# Install RDKit from conda-forge channel
RUN conda install conda-forge::rdkit
# Setting working directory
WORKDIR /src
# Copying requirements file
COPY requirements.txt requirements.txt
# Python dependencies installation while reducing size of an image
RUN pip install --no-cache-dir -r requirements.txt
# Port exposion

EXPOSE 8000
# Entry point for container
ENTRYPOINT [ "uvicorn" ]
# Run app
CMD [ "main:app", "--host","0.0.0.0","--port","8000" ]