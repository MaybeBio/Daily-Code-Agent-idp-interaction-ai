FROM vemlp-cn-beijing.cr.volces.com/preset-images/pytorch:2.7.1-cu12.6.3-py3.11-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Asia/Shanghai \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CUTLASS_PATH=/opt/cutlass

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    g++ \
    gcc \
    libc6-dev \
    make \
    postgresql \
    hmmer \
    kalign \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Clone CUTLASS
RUN git clone -b v3.5.1 https://github.com/NVIDIA/cutlass.git /opt/cutlass

# Clone Protenix and install it (this also pulls its requirements.txt)
RUN git clone https://github.com/bytedance/Protenix.git /app \
    && pip3 install --no-cache-dir -r /app/requirements.txt -i https://pypi.org/simple \
    && pip3 install --no-cache-dir -e /app

# Set default entrypoint
ENTRYPOINT ["python", "/app/runner/inference.py"]
