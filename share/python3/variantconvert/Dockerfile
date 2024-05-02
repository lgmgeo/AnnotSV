FROM centos:7
LABEL org.opencontainers.image.authors="https://github.com/SamuelNicaise/variantconvert"
LABEL Software="variantconvert"

RUN yum -y install sudo make wget bzip2 gcc git

#conda
RUN cd /tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/
RUN echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh
RUN chmod 755 /etc/profile.d/miniconda.sh
ENV PATH="$PATH:/usr/local/lib/miniconda3/bin"
RUN chmod -R 755 /usr/local/lib/miniconda3/
RUN rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh

#work env
RUN conda create -n common python=3.10 pandas

#bcftools
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.14
RUN yum -y install zlib-devel xz-devel bzip2-devel libcurl-devel
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2
WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install
WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION

#htslib
ENV HTSLIB_INSTALL_DIR=/opt/htslib
ENV HTSLIB_VERSION=1.14
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
  tar --bzip2 -xf htslib-$HTSLIB_VERSION.tar.bz2
WORKDIR /tmp/htslib-$HTSLIB_VERSION
RUN make prefix=$HTSLIB_INSTALL_DIR && \
  make prefix=$HTSLIB_INSTALL_DIR install
WORKDIR /
RUN ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip
RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix
RUN rm -rf /tmp/htslib-$HTSLIB_VERSION

#install varianconvert
RUN git clone https://github.com/SamuelNicaise/variantconvert.git 
WORKDIR variantconvert/
RUN /usr/local/lib/miniconda3/envs/common/bin/pip install -e .

RUN echo "source activate common" > ~/.bashrc
ENTRYPOINT ["/usr/local/lib/miniconda3/envs/common/bin/python", "/variantconvert/variantconvert/__main__.py"]
