FROM ubuntu:16.04

RUN apt-get clean all && apt-get update && apt-get install -y build-essential apt-utils git wget perl \
    python3.5 python2.7 python3-pip python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2 exonerate \
    bzip2 liblzma-dev libbz2-dev

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN echo "mysql-server mysql-server/root_password password lorean" | debconf-set-selections

RUN echo "mysql-server mysql-server/root_password_again password lorean" | debconf-set-selections

RUN apt-get install -y mysql-server mysql-client mysql-common bowtie bioperl apache2 libcairo2-dev libpango1.0-dev 

RUN pip3 install biopython==1.68 bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils regex pysam matplotlib progressbar2 \
    psutil memory_profiler pathlib colorama

WORKDIR /opt/

RUN git clone git://github.com/pezmaster31/bamtools.git && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  sudo ln -f -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .

RUN git clone -b LoReAn_ipscan --single-branch https://github.com/lfaino/LoReAn.git  && git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git && \
    mv Complete-Striped-Smith-Waterman-Library SW && cd SW/src/ && make && cp ssw_lib.py  /opt/LoReAn/code/ && cp libssw.so  /opt/LoReAn/code/

WORKDIR /opt/LoReAn/third_party/software/

RUN tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean &&\
 sudo ./configure --prefix=$PWD && sudo make && sudo make install

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz && tar -zxvf gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz

RUN wget https://github.com/PASApipeline/PASApipeline/archive/v2.1.0.tar.gz && tar -zxvf v2.1.0.tar.gz && rm v2.1.0.tar.gz &&\
    mv PASApipeline-2.1.0 PASApipeline && cd PASApipeline && make clean && make && cd .. &&  cp ../conf_files/conf.txt PASApipeline/pasa_conf/ &&\
    cp ../scripts/process_GMAP_alignments_gff3_chimeras_ok.pl PASApipeline/scripts/ &&\
    chmod 775 PASApipeline/scripts/process_GMAP_alignments_gff3_chimeras_ok.pl

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz && \
    tar -zxvf augustus.current.tar.gz && rm augustus.current.tar.gz && cd augustus  && make clean && make

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.5.1.tar.gz && tar -zxvf Trinity-v2.5.1.tar.gz && \
    mv trinityrnaseq-Trinity-v2.5.1 Trinity &&rm Trinity-v2.5.1.tar.gz && cd Trinity && make && make plugins

RUN git clone https://github.com/alexdobin/STAR.git

RUN wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz &&  tar -zxvf v3.0.1.tar.gz && rm v3.0.1.tar.gz &&\
    cd TransDecoder-3.0.1 && make

RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-06-20.tar.gz && tar -zxvf gmap-gsnap-2017-06-20.tar.gz && rm gmap-gsnap-2017-06-20.tar.gz && \
    mv gmap-2017-06-20/ gmap && cd gmap/ && ./configure && make && sudo make install

RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8e.tar.gz && tar -zxvf fasta-36.3.8e.tar.gz && rm fasta-36.3.8e.tar.gz &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 && cp /opt/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /usr/local/bin/fasta

RUN wget http://exon.gatech.edu/Braker/BRAKER2.tar.gz && tar -zxvf BRAKER2.tar.gz && rm BRAKER2.tar.gz && mv BRAKER* BRAKER && \
 wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i  Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric && sudo cpan -f -i File::Which

RUN mkdir gffread && cd gffread && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /usr/local/bin

RUN wget http://genometools.org/pub/genometools-1.5.9.tar.gz && tar -zxvf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make

RUN cp ../../code/createUser.py /usr/local/bin && chmod 775 /usr/local/bin/createUser.py

RUN cp ../conf_files/pathToExport.txt /etc/profile.d/pathToExport.sh

RUN cp ../conf_files/extrinsic.M.RM.E.W.P.cfg augustus/config/extrinsic/

RUN rm /opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl

RUN sudo chmod -R 775 /opt/LoReAn/code/

RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz.md5 && \
    md5sum -c interproscan-5.27-66.0-64-bit.tar.gz.md5

RUN tar -pxvzf interproscan-5.27-66.0-64-bit.tar.gz && rm interproscan-5.27-66.0-64-bit.tar.gz

WORKDIR /opt/LoReAn/third_party/software/interproscan-5.27-66.0/data

RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz.md5
    md5sum -c panther-data-12.0.tar.gz.md5

RUN tar -pxvzf panther-data-12.0.tar.gz && rm panther-data-12.0.tar.gz.md5

WORKDIR /opt/LoReAn/third_party/software/interproscan-5.27-66.0/

#COPY signalp-4.1f.Linux.tar.gz /
#RUN  tar -xzf /signalp-4.1f.Linux.tar.gz -C /interproscan-5.27-66.0/bin/signalp/4.1 --strip-components 1

#COPY tmhmm-2.0c.Linux.tar.gz /
#RUN  tar -xzf /tmhmm-2.0c.Linux.tar.gz -C / && \
#     cp /tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64  /interproscan-5.27-66.0/bin/tmhmm/2.0c/decodeanhmm && \
#     cp /tmhmm-2.0c/lib/TMHMM2.0.model  /interproscan-5.27-66.0/data/tmhmm/2.0c/TMHMM2.0c.model

WORKDIR /data/

CMD /usr/local/bin/createUser.py