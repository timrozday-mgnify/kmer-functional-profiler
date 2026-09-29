# kMermaid (no conda or PyPI package) at a fixed commit, with the benchmark's wrapper
# (kmermaid_kfp.py: retrain on KO clusters, classify reads). Built by containers/build.sh.
FROM python:3.11-slim
ARG KMERMAID_COMMIT=edcb4ed030ff16b0c5005e3d704a9da0fb61e049
RUN apt-get update && apt-get install -y --no-install-recommends git procps \
    && pip install --no-cache-dir numpy pandas pyarrow \
        "git+https://github.com/AuslanderLab/kmermaid@${KMERMAID_COMMIT}" \
    && apt-get purge -y git && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*
COPY kmermaid_kfp.py /usr/local/bin/kmermaid_kfp.py
RUN sed -i '1i #!/usr/bin/env python3' /usr/local/bin/kmermaid_kfp.py \
    && chmod +x /usr/local/bin/kmermaid_kfp.py
