FROM us-docker.pkg.dev/deeplearning-platform-release/gcr.io/tf-cu113.2-10.py310:latest

RUN apt-get update && apt-get install -y git

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
        git \
        time \
        make \
        gcc \
        g++ \
        curl \
        bzip2 \
        build-essential \
        libdb-dev \
        liburi-perl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


RUN python3 -m pip install --upgrade \
    pip setuptools wheel \
    "hatchling>=1.26" \
    "packaging>=24.0"

RUN pip install "pyfamsa<0.6.0"

RUN pip install pyBigWig bio scikit-learn biopython bcbio-gff requests


RUN cd /opt && \
    git        clone https://github.com/Gaius-Augustus/Tiberius && \
    cd      Tiberius && \
    pip install . && \
    chmod +x tiberius.py && \
    chmod +x tiberius/*py

RUN mkdir -p /opt/Tiberius/model_weights && chmod -R 777 /opt/Tiberius/model_weights

ENV PATH=${PATH}:/opt/Tiberius/tiberius/
ENV PATH=${PATH}:/opt/Tiberius/


# micromamba installation
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C /usr/local/bin --strip-components=1 bin/micromamba

RUN micromamba create -y -n bio -c conda-forge -c bioconda \
        minimap2 \
        stringtie \
        samtools \
        transdecoder \
        diamond \
        bedtools \
        augustus \
    && micromamba clean -a -y

RUN  cd /opt && \
     git clone https://github.com/TransDecoder/TransDecoder
ENV PATH=${PATH}:/opt/TransDecoder/util

RUN  cd /opt && \
     git clone https://github.com/tomasbruna/miniprothint
ENV PATH=${PATH}:/opt/miniprothint

RUN cd /opt && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && make
ENV PATH=${PATH}:/opt/miniprot

RUN  cd /opt && \
     git clone https://github.com/tomasbruna/miniprot-boundary-scorer && \
     cd miniprot-boundary-scorer && make
ENV PATH=${PATH}:/opt/miniprot-boundary-scorer


RUN cd /opt && \
        wget -O hisat2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
        unzip hisat2.zip

ENV PATH=${PATH}:/opt/hisat2-2.2.1/

RUN ln -s /opt/micromamba/envs/bio/bin/minimap2     /usr/local/bin/minimap2     && \
    ln -s /opt/micromamba/envs/bio/bin/stringtie    /usr/local/bin/stringtie    && \
    ln -s /opt/micromamba/envs/bio/bin/samtools     /usr/local/bin/samtools     && \
    ln -s /opt/micromamba/envs/bio/bin/TransDecoder.LongOrfs  /usr/local/bin/TransDecoder.LongOrfs  && \
    ln -s /opt/micromamba/envs/bio/bin/TransDecoder.Predict   /usr/local/bin/TransDecoder.Predict   && \
    ln -s /opt/micromamba/envs/bio/bin/diamond      /usr/local/bin/diamond      && \
    ln -s /opt/micromamba/envs/bio/bin/bedtools     /usr/local/bin/bedtools     && \
    ln -s /opt/micromamba/envs/bio/bin/bam2hints    /usr/local/bin/bam2hints


RUN micromamba install -y -n bio -c bioconda -c conda-forge sra-tools \
    && micromamba clean -a -y
RUN ln -s /opt/micromamba/envs/bio/bin/fasterq-dump /usr/local/bin/fasterq-dump
RUN ln -s /opt/micromamba/envs/bio/bin/prefetch /usr/local/bin/prefetch
RUN ln -s /opt/micromamba/envs/bio/bin/vdb-config /usr/local/bin/vdb-config


RUN micromamba install -c bioconda -c conda-forge seqkit

ENTRYPOINT ["bash"]


USER ${NB_UID}