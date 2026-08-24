FROM ubuntu:20.04

LABEL maintainer="Yoshitaka Sakamoto"

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    curl \
    ca-certificates \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH $PATH:/root/.cargo/bin

WORKDIR /opt

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar jxvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17/htslib-1.17 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

# Install BWA
RUN wget -O bwa-0.7.18.tar.gz https://github.com/lh3/bwa/archive/refs/tags/v0.7.18.tar.gz && \
    tar xzf bwa-0.7.18.tar.gz && \
    cd bwa-0.7.18 && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa-0.7.18*

# Install nanomonsv_postprocess. VERSION is the git tag (or branch) to build, so
# no version string has to be edited here per release; declared last so every
# layer above stays cached when it changes.
ARG VERSION
RUN test -n "${VERSION}" || { \
        echo "VERSION build-arg is required, e.g. --build-arg VERSION=v0.3.0" >&2; \
        exit 1; \
    }
RUN git clone --depth 1 --branch "${VERSION}" \
        https://github.com/yos-sk/nanomonsv_postprocess.git \
        /opt/nanomonsv_postprocess

WORKDIR /opt/nanomonsv_postprocess
# Tags that carry a Cargo.lock build exactly what was tested; older tags predate
# it and need the hts-sys pin resolved here instead (a fresh resolve picks up
# hts-sys 2.2.1, whose size_t bindings no longer match rust-htslib 0.46).
RUN if [ -f Cargo.lock ]; then \
        cargo build --release --locked; \
    else \
        cargo generate-lockfile && \
        cargo update -p hts-sys --precise 2.1.4 && \
        cargo build --release; \
    fi

ENV PATH="/opt/nanomonsv_postprocess/target/release:$PATH"

# Fail the build rather than push a broken image.
RUN nanomonsv_postprocess --version && bwa 2>&1 | head -1 && \
    samtools --version | head -1

CMD ["/bin/bash"]
