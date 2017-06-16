
FROM google/cloud-sdk

LABEL maintainer="R. Jay Mashl <rmashl@wustl.edu>"
LABEL program="PanCanAtlas analysis"
LABEL version="0.1"

RUN apt-get update && apt-get -y install \
    autoconf \
    build-essential \
    libncurses-dev \
    perl \    		   
    pkg-config \
    unzip \
    wget \
    zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*
		 
# install vcfanno
WORKDIR /usr/local/bin
RUN    wget -O vcfanno https://github.com/brentp/vcfanno/releases/download/v0.2.6/vcfanno_linux64
RUN    chmod +x ./vcfanno

# install vcftools
WORKDIR /usr/local/src
RUN    wget -O v0.1.14.zip https://github.com/vcftools/vcftools/archive/v0.1.14.zip && unzip v0.1.14.zip && rm -f v0.1.14.zip  && cd vcftools-0.1.14 && export ZLIB_LIBS=-lz && export ZLIB_CFLAGS=-I/usr/include && ./autogen.sh && ./configure --prefix=/usr/local && make && make install

# install samtools
WORKDIR /usr/local/src
RUN     wget -O samtools-1.2.tar.bz2 https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && tar xjf samtools-1.2.tar.bz2  && rm -f samtools-1.2.tar.bz2
RUN     cd samtools-1.2/htslib-1.2.1 && ./configure && make  && cp bgzip htsfile tabix /usr/local/bin/ && cp libhts.so.1 /usr/local/lib/ && /sbin/ldconfig
WORKDIR /usr/local/src
RUN     cd samtools-1.2 && make  && cp samtools /usr/local/bin/

# scripts
WORKDIR /usr/local/bin
COPY     variant_QC_annotation.sh  filter_VCF_AF_AD.py ExAC_config.toml  ./
RUN      chmod +x ./variant_QC_annotation.sh


