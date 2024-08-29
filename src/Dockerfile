FROM python:3.9-slim
# Use Minicoda image  
FROM continuumio/miniconda3
# Install RDKit from conda-forge channel
RUN conda install conda-forge::rdkit
# Install psycog2-binary
RUN pip install psycopg2-binary

# Setting working directory
WORKDIR /app


# Copying requirements file
COPY requirements.txt requirements.txt
# Python dependencies installation while reducing size of an image
RUN pip install --no-cache-dir -r requirements.txt
# Port exposion

EXPOSE 8000
# Entry point for container
ENTRYPOINT [ "uvicorn" ]

ENV PYTHONPATH=/src
# Run app
CMD [ "main:app", "--host","0.0.0.0","--port","8000" ]