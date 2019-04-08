FROM ubuntu:16.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y
RUN apt-get install -y --no-install-recommends \
        apt-transport-https \
        build-essential \
        git \
        libcurl4-openssl-dev \
        libxml2-dev \
        software-properties-common \
        libreadline-dev \
        libpcre++-dev \
        libblas-dev \
        liblapack-dev \
        libatlas-base-dev \
        gfortran \
        liblzma-dev \
        libbz2-dev \
        locales \
        zlib1g-dev \
        xserver-xorg \ 
        xdm \
        xfonts-base \ 
        xfonts-100dpi \
        xfonts-75dpi \
        libqt5x11extras5

# Locales
RUN locale-gen en_GB.UTF-8
ENV LANG en_GB.UTF-8
ENV LANGUAGE en_GB:en
ENV LC_ALL en_GB.UTF-8

# Python 3.6
RUN add-apt-repository -y ppa:deadsnakes/ppa && \
        apt remove -y python3-apt && \
        apt-get update && \
        apt-get install -y --no-install-recommends \
                python3.6 python3.6-dev python3-pip python3-apt python3-tk \ 
                python3-setuptools \
                software-properties-common python3-software-properties && \
        rm /usr/bin/python3 && \
        ln -s /usr/bin/python3.6 /usr/bin/python3 && \
        cp /usr/lib/python3/dist-packages/apt_pkg.cpython-35m-x86_64-linux-gnu.so \
            /usr/lib/python3/dist-packages/apt_pkg.so

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.6 1
RUN update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

# R 3.3
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#    add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/' && \
#    apt-get update && \
#    apt-get install -y --no-install-recommends \
#        r-recommended=3.3.3-1xenial0 r-base=3.3.3-1xenial0

#RUN R -e "install.packages(c('robustbase', 'rrcov'), repos='http://cran.us.r-project.org')"

#RUN R -e "source('http://bioconductor.org/biocLite.R'); \
#         biocLite(''); \
#         biocLite('GEOquery'); \
#         biocLite('data.table'); \
#         biocLite('rhdf5'); \
#         biocLite('limma'); \
#         biocLite('oligo'); \
#         biocLite('ArrayExpress'); \
#         biocLite('magrittr'); \
#         biocLite('dplyr'); \
#         biocLite('AnnotationDbi'); \
#         biocLite('WGCNA');"

# RUN R -e "install.packages(c('INSPIRE', 'MASS', 'lme4', 'lmerTest'), repos='http://cran.us.r-project.org')"

RUN pip3 install --upgrade pip
RUN pip3 install wheel
RUN pip3 install numpy scipy docutils six pytest matplotlib lxml PyQt5 ete3

COPY . /app
WORKDIR /app
RUN python3 setup.py install
