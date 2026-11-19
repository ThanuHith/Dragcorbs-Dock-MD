FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    HOME=/app \
    PORT=8080

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates libxrender1 libxext6 libexpat1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --upgrade pip && pip install \
    fastapi "uvicorn[standard]" pydantic numpy pandas requests \
    rdkit meeko gemmi streamlit

COPY backend /app/backend
COPY frontend /app/frontend
COPY start.sh /app/start.sh

RUN rm -f /app/backend/vina.exe && \
    curl -fsSL "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64" \
    -o /app/backend/vina && \
    chmod +x /app/backend/vina /app/start.sh && \
    chmod -R 777 /app && \
    /app/backend/vina --version

EXPOSE 8080

CMD ["/app/start.sh"]
