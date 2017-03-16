FROM ubuntu:16.04

RUN apt-get clean all && apt-get update && apt-get install -y build-essential apt-utils git wget perl \
    python3.5 python-pip3 debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN echo "mysql-server mysql-server/root_password password lorean" | debconf-set-selections
RUN echo "mysql-server mysql-server/root_password_again password lorean" | debconf-set-selections

RUN apt-get install -y mysql-server mysql-client mysql-common bowtie bioperl apache2 libcairo2-dev libpango1.0-dev 

RUN pip3 install biopython==1.68 bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils
          
RUN adduser --disabled-password --gecos '' lorean &&\
    adduser lorean sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER lorean

WORKDIR /home/lorean

RUN mkdir bin

WORKDIR /home/lorean/bin

RUN git clone git://github.com/pezmaster31/bamtools.git && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  sudo ln -f -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .

RUN git clone https://github.com/lfaino/LoReAn.git  

WORKDIR /home/lorean/bin/LoReAn/third_party/software/

RUN tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean && sudo ./configure --prefix=$PWD && sudo make && sudo make install 

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz

RUN tar -zxvf gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz

RUN wget https://github.com/PASApipeline/PASApipeline/archive/v2.0.2.tar.gz && tar -zxvf v2.0.2.tar.gz && rm v2.0.2.tar.gz &&\
    cd PASApipeline-2.0.2 && make clean && make && cd .. &&  cp ../conf_files/conf.txt PASApipeline-2.0.2/pasa_conf/ 
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz && \
    tar -zxvf augustus-3.2.3.tar.gz && rm augustus-3.2.3.tar.gz && mv augustus-3.2.3 augustus && cd augustus  && make clean && make 
    
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz && tar -zxvf Trinity-v2.4.0.tar.gz && mv trinityrnaseq-Trinity-v2.4.0 Trinity && \
    rm Trinity-v2.4.0.tar.gz && cd Trinity && make && make plugins 

RUN wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz && tar -xzf 2.5.2b.tar.gz && rm 2.5.2b.tar.gz &&\
    cd STAR-2.5.2b/source && make STAR && git submodule update --init --recursive 

RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2015-12-31.v10.tar.gz && tar -zxvf gmap-gsnap-2015-12-31.v10.tar.gz && rm gmap-gsnap-2015-12-31.v10.tar.gz &&\
    mv gmap-2015-12-31// gmap && cd gmap/  && ./configure && make && sudo make install 

RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8e.tar.gz && tar -zxvf fasta-36.3.8e.tar.gz && rm fasta-36.3.8e.tar.gz &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 && cp /home/lorean/bin/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /home/lorean/bin/fasta

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/BRAKER1.tar.gz && tar -zxvf BRAKER1.tar.gz && rm BRAKER1.tar.gz

RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz 
   
   
RUN sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i  Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric 
     
RUN mkdir gffread && cd gffread && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /home/lorean/bin

RUN wget http://genometools.org/pub/genometools-1.5.9.tar.gz && \
     tar -zxvf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make
     
RUN cat ~/.bashrc ../conf_files/pathToExport.txt > ~/.bashrc_new && mv ~/.bashrc_new ~/.bashrc && source ~/.bashrc 

    
WORKDIR /data/
