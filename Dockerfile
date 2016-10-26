FROM ubuntu:16.04

RUN apt-get clean all && apt-get update && apt-get install -y build-essential apt-utils git wget perl \
    python2.7 python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN echo "mysql-server mysql-server/root_password password lorean" | debconf-set-selections
RUN echo "mysql-server mysql-server/root_password_again password lorean" | debconf-set-selections

RUN apt-get install -y mysql-server mysql-client mysql-common bowtie bioperl apache2 libcairo2-dev libpango1.0-dev 

RUN pip install biopython bcbio-gff
          
RUN adduser --disabled-password --gecos '' lorean &&\
    adduser lorean sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER lorean

WORKDIR /home/lorean

RUN mkdir bin

WORKDIR /home/lorean/bin

RUN git clone git://github.com/pezmaster31/bamtools.git && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  ln -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  ln -s /usr/local/lib/bamtools/libbamtools.* .

RUN git clone https://github.com/lfaino/LoReAn.git  

WORKDIR /home/lorean/bin/LoReAn/third_party/software/

RUN cd AATpackage-r03052011 && make clean && sudo ./configure --prefix=$PWD && sudo make && sudo make install &&\
    cd .. && cd /home/lorean/bin/LoReAn/third_party/software/ && cd PASApipeline-2.0.2 && make clean && make && cd .. &&\
    cp ../conf_files/conf.txt PASApipeline-2.0.2/pasa_conf/ && cd /home/lorean/bin/LoReAn/third_party/software/augustus && make clean && make &&\
    cd /home/lorean/bin/LoReAn/third_party/software/ && cd trinityrnaseq-2.2.0 && make && make plugins && cd /home/lorean/bin/LoReAn/third_party/software/ &&\
    cd STAR-2.5.0b && make STAR && git submodule update --init --recursive && cd /home/lorean/bin/LoReAn/third_party/software/  && \
    cd gmap-2015-12-31/ && ./configure && make && make install && cd /home/lorean/bin/LoReAn/third_party/software/ &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 &&\
    cp /home/lorean/bin/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /home/lorean/bin/fasta &&\
    cd /home/lorean/bin/LoReAn/third_party/software/  && cd genometools-1.5.9 && make &&\
    cd /home/lorean/bin/LoReAn/third_party/software/ &&\
    cat ~/.bashrc ../conf_files/pathToExport.txt > ~/.bashrc_new && mv ~/.bashrc_new ~/.bashrc && source ~/.bashrc && \
    cp ../conf_files/gm_key ~/.gm_key 
    
RUN sudo perl -MCPAN -e shell && cpan -f -i YAML && cpan -f -i Hash::Merge && cpan -f -i  Logger::Simple && cpan -f -i  Parallel::ForkManager &&\
    cpan -f -i Config::Std && cpan -f -i Scalar::Util::Numeric 
     
RUN mkdir gffreadF && cd gffreadF && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /home/lorean/bin
   
WORKDIR /data/
