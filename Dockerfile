# LunaSRIM web application
#   docker build -t lunasrim .
#   docker run -p 8000:8000 lunasrim        ->  http://localhost:8000/
FROM python:3.12-slim
WORKDIR /app
COPY . /app
RUN sh tools/fetch_web_vendor.sh || echo "plotly will be loaded from the CDN" \
 && pip install --no-cache-dir ".[web]"
ENV LUNASRIM_MAX_IONS=100000 LUNASRIM_MAX_JOBS=2
EXPOSE 8000
CMD ["python", "-m", "lunasrim.web", "--host", "0.0.0.0", "--port", "8000"]
