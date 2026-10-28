FROM python:3.11-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    gcc \
    g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Clone and install psauron
RUN git clone  https://github.com/salzberg-lab/PSAURON.git /app \
    && pip install --no-cache-dir .

# Set entrypoint to the psauron CLI
ENTRYPOINT ["psauron"]
CMD ["--help"]
