# HUMAnN 3.9 biocontainer (MetaPhlAn 4.1.1, DIAMOND 2.0.15) with bowtie2 2.5.5: the
# biocontainer resolves bowtie2 2.2.3, whose bowtie2-build has no --threads, so HUMAnN's
# index build fails. Built by containers/build.sh; the HUMAnN 4 image builds on it.
FROM quay.io/biocontainers/humann:3.9--py312hdfd78af_0
ARG BOWTIE2=2.5.5
RUN python -c "import io, urllib.request, zipfile; zipfile.ZipFile(io.BytesIO(urllib.request.urlopen('https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2}/bowtie2-${BOWTIE2}-linux-x86_64.zip').read())).extractall('/opt')" \
    && chmod +x /opt/bowtie2-${BOWTIE2}-linux-x86_64/bowtie2* \
    && ln -sf /opt/bowtie2-${BOWTIE2}-linux-x86_64/bowtie2* /usr/local/bin/
