FROM python:3.12-slim

# Copy pyproject.toml (and source code if needed)
COPY pyproject.toml /app/pyproject.toml
COPY . /app

WORKDIR /app

# Install directly from pyproject.toml
RUN python -m venv /opt/venv && \
    /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install .

# Environment
ENV PATH=/opt/venv/bin:$PATH
ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["/bin/sh", "-c"]