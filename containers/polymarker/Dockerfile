FROM almalinux:9

ARG POLYMARKER_VERSION=1.3.3
ARG MAFFT_VERSION=7.526
ARG BLAST_VERSION=2.16.0
ARG PRIMER3_VERSION=2.6.1
ARG EXONERATE_VERSION=2.4.0

LABEL org.opencontainers.image.source="https://github.com/JohGi/KASPberry"
LABEL org.opencontainers.image.title="KASPberry PolyMarker runtime"
LABEL org.opencontainers.image.description="Reproducible PolyMarker runtime used by KASPberry for KASP assay design."
LABEL org.opencontainers.image.version="${POLYMARKER_VERSION}"
LABEL org.opencontainers.image.authors="Johanna Girodolle"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Build/runtime dependencies required by PolyMarker and its external tools.
RUN dnf -y group install "Development Tools" && \
    dnf -y install \
        wget \
        ruby \
        glib2-devel \
        ncurses-devel \
        bzip2-devel && \
    dnf clean all && \
    rm -rf /var/cache/dnf

WORKDIR /tmp/build

# MAFFT
RUN wget -q \
        "https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-gcc_fc6.x86_64.rpm" && \
    dnf -y install \
        "./mafft-${MAFFT_VERSION}-gcc_fc6.x86_64.rpm" && \
    rm -f "mafft-${MAFFT_VERSION}-gcc_fc6.x86_64.rpm" && \
    dnf clean all && \
    rm -rf /var/cache/dnf

# NCBI BLAST+
RUN wget -q \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-1.x86_64.rpm" && \
    dnf -y install \
        "./ncbi-blast-${BLAST_VERSION}+-1.x86_64.rpm" && \
    rm -f "ncbi-blast-${BLAST_VERSION}+-1.x86_64.rpm" && \
    dnf clean all && \
    rm -rf /var/cache/dnf

# Primer3
RUN wget -q \
        "https://github.com/primer3-org/primer3/archive/refs/tags/v${PRIMER3_VERSION}.tar.gz" \
        -O "primer3-${PRIMER3_VERSION}.tar.gz" && \
    tar -xzf "primer3-${PRIMER3_VERSION}.tar.gz" && \
    make -C "primer3-${PRIMER3_VERSION}/src" && \
    make -C "primer3-${PRIMER3_VERSION}/src" install && \
    rm -rf \
        "primer3-${PRIMER3_VERSION}" \
        "primer3-${PRIMER3_VERSION}.tar.gz"

# Exonerate
RUN wget -q \
        "https://github.com/cb2e6f/exonerate/archive/refs/tags/v${EXONERATE_VERSION}.tar.gz" \
        -O "exonerate-${EXONERATE_VERSION}.tar.gz" && \
    tar -xzf "exonerate-${EXONERATE_VERSION}.tar.gz" && \
    cd "exonerate-${EXONERATE_VERSION}" && \
    ./configure && \
    make && \
    make install && \
    cd /tmp/build && \
    rm -rf \
        "exonerate-${EXONERATE_VERSION}" \
        "exonerate-${EXONERATE_VERSION}.tar.gz"

# PolyMarker.
#
# bio-polymarker 1.3.3 uses SortedSet at runtime. With Ruby 3, SortedSet is
# provided by the separate sorted_set gem, whose rbtree dependency contains
# a native extension and therefore requires the Ruby development headers.
ARG SORTED_SET_VERSION=1.0.3

RUN dnf -y install ruby-devel && \
    dnf clean all && \
    rm -rf /var/cache/dnf

RUN gem install --no-document rake && \
    gem install --no-document sorted_set \
        --version "${SORTED_SET_VERSION}" && \
    gem install --no-document bio-polymarker \
        --version "${POLYMARKER_VERSION}"

# Build-time smoke tests.
# Keep each dependency check separate so failures are easy to diagnose.

RUN ruby --version && \
    gem list --local bio-polymarker \
        | grep -F "bio-polymarker (${POLYMARKER_VERSION})" && \
    gem list --local sorted_set \
        | grep -F "sorted_set (${SORTED_SET_VERSION})" && \
    ruby -e 'require "set"; require "sorted_set"; s = SortedSet.new([3,1,2]); abort unless s.to_a == [1,2,3]' && \
    command -v polymarker.rb >/dev/null && \
    polymarker.rb --help >/dev/null

RUN mafft --version 2>&1

RUN primer3_core --about 2>&1

# Exonerate prints the requested version correctly but returns exit code 1
# for --version. Validate the output explicitly instead.
RUN exonerate_version="$(exonerate --version 2>&1 || true)" && \
    echo "${exonerate_version}" && \
    grep -Fq "exonerate version ${EXONERATE_VERSION}" <<< "${exonerate_version}"

RUN blastn -version

WORKDIR /

RUN rm -rf /tmp/build
