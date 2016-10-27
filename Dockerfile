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

WORKDIR /home/lorean/bin/LoReAn/third_party/

RUN tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean && sudo ./configure --prefix=$PWD && sudo make && sudo make install 

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz

RUN wget https://github.com/PASApipeline/PASApipeline/archive/v2.0.2.tar.gz && tar -zxvf v2.0.2.tar.gz && rm v2.0.2.tar.gz &&\
    cd PASApipeline-2.0.2 && make clean && make && cd .. &&  cp ../conf_files/conf.txt PASApipeline-2.0.2/pasa_conf/ 
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.2.tar.gz && mkdir augustus && cd augustus && mv ../augustus-3.2.2.tar.gz . &&\
    tar -zxvf augustus-3.2.2.tar.gz && rm augustus-3.2.2.tar.gz  && make clean && make 
    
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.2.0.tar.gz && tar -zxvf v2.2.0.tar.gz && rm v2.2.0.tar.gz && cd trinityrnaseq-2.2.0 && make && make plugins 

RUN wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz && tar -xzf 2.5.2b.tar.gz && rm 2.5.2b.tar.gz &&\
    cd STAR-2.5.2b && make STAR && git submodule update --init --recursive 

RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2016-09-23.tar.gz && tar -zxvf gmap-gsnap-2016-09-23.tar.gz && rm gmap-gsnap-2016-09-23.tar.gz &&\
    mv gmap-2016-09-23/ gmap && cd gmap/ && ./configure && make && make install 

RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8e.tar.gz && tar -zxvf fasta-36.3.8e.tar.gz && rm fasta-36.3.8e.tar.gz &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 && cp /home/lorean/bin/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /home/lorean/bin/fasta

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/BRAKER1.tar.gz && tar -zxvf BRAKER1.tar.gz

RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz 
   
   
RUN sudo perl -MCPAN -e shell && cpan -f -i YAML && cpan -f -i Hash::Merge && cpan -f -i  Logger::Simple && cpan -f -i  Parallel::ForkManager &&\
    cpan -f -i Config::Std && cpan -f -i Scalar::Util::Numeric 
     
RUN mkdir gffreadF && cd gffreadF && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /home/lorean/bin

RUN wget http://genometools.org/pub/genometools-1.5.9.tar.gz && \
     tar -zxvf genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make
RUN cat ~/.bashrc ../conf_files/pathToExport.txt > ~/.bashrc_new && mv ~/.bashrc_new ~/.bashrc && source ~/.bashrc && \
    cp ../conf_files/gm_key ~/.gm_key
    
WORKDIR /data/
