#!/usr/bin/env bash
set -e

echo "Starting Dragcorbs FastAPI Physics Backend on port 8000..."
cd /app/backend
uvicorn main:app --host 127.0.0.1 --port 8000 &

for i in $(seq 1 30); do
  if curl -s http://127.0.0.1:8000/health > /dev/null; then
    echo "Dragcorbs Backend is healthy!"
    break
  fi
  sleep 1
done

echo "Starting Dragcorbs Streamlit Web App on port ${PORT:-8080}..."
cd /app/frontend
export DRAGCORBS_BACKEND_URL="http://127.0.0.1:8000"
exec streamlit run app.py \
  --server.port "${PORT:-8080}" \
  --server.address "0.0.0.0" \
  --server.headless true \
  --browser.gatherUsageStats false