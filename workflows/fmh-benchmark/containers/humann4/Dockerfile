# HUMAnN 4 (4.0.0.alpha.2, GitHub only: not on bioconda or PyPI) at a fixed commit, over the
# HUMAnN 3.9 image (containers/humann3), which has the MetaPhlAn (4.1.1) and DIAMOND
# (2.0.15) versions HUMAnN 4 pins, and bowtie2 2.5.5. Built by containers/build.sh.
FROM kfp-humann3:3.9-bt2.5.5
ARG HUMANN_COMMIT=e07b3a34d0b94c09a8ac5d28ff95009611178be2
RUN (python -m pip --version || python -m ensurepip) \
    && python -m pip install --no-cache-dir --no-deps \
        "https://github.com/biobakery/humann/archive/${HUMANN_COMMIT}.tar.gz"
