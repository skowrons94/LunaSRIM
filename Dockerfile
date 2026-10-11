# LunaSRIM web application
#   docker build -t lunasrim .
#   docker run -p 8000:8000 lunasrim        ->  http://localhost:8000/
FROM python:3.12-slim
WORKDIR /src
COPY pyproject.toml MANIFEST.in README.md LICENSE CHANGELOG.md CITATION.cff ./
COPY lunasrim ./lunasrim
COPY tools/fetch_web_vendor.sh ./tools/
# local copies of plotly.js, marked and KaTeX (otherwise the browser loads them from a CDN)
RUN (sh tools/fetch_web_vendor.sh || echo "the JavaScript libraries will be loaded from a CDN") \
 && pip install --no-cache-dir ".[web]" \
 && cd / && rm -rf /src \
 && useradd --create-home lunasrim
USER lunasrim
WORKDIR /home/lunasrim
ENV LUNASRIM_MAX_IONS=100000 LUNASRIM_MAX_JOBS=2
EXPOSE 8000
CMD ["lunasrim", "web", "--host", "0.0.0.0", "--port", "8000"]
