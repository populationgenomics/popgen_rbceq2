# Use a slim Python image for efficiency
FROM python:3.12-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y procps git && rm -rf /var/lib/apt/lists/*

## (Nextflow's bin/ script is automatically added to the PATH by Nextflow,
#  so we don't strictly need to COPY it into the container if we mount it,
#  but for a truly self-contained container on GCP AR, we can include it.)
RUN pip install --no-cache-dir "git+https://github.com/limcintyre/RBCeq2.git@561d8a3362f45a86a543f453b1251cd8998bb743"

# Standard entrypoint (though Nextflow overrides this)
CMD ["python3"]
