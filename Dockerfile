FROM ubuntu:xenial

RUN apt-get clean all && apt-get update && apt-get install -y -q build-essential git wget perl \
    python3.5 python2.7 software-properties-common python3-pip python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2 exonerate \
    bzip2 liblzma-dev libbz2-dev software-properties-common libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev cmake tcsh g++ iputils-ping


RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN echo "mysql-server mysql-server/root_password password lorean" | debconf-set-selections

RUN echo "mysql-server mysql-server/root_password_again password lorean" | debconf-set-selections

RUN apt-get install -y mysql-server mysql-client mysql-common bowtie bioperl apache2 libcairo2-dev libpango1.0-dev 

RUN pip3 install numpy biopython==1.68 bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils regex pysam matplotlib progressbar2 \
    psutil memory_profiler pathlib colorama simplesam

WORKDIR /opt/

RUN git clone git://github.com/pezmaster31/bamtools.git && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  sudo ln -f -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .

RUN git clone -b dev --single-branch https://github.com/lfaino/LoReAn.git  && git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git && \
    mv Complete-Striped-Smith-Waterman-Library SW && cd SW/src/ && make && cp ssw_lib.py  /opt/LoReAn/code/ && cp libssw.so  /opt/LoReAn/code/

WORKDIR /opt/LoReAn/third_party/software/

RUN tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean &&\
 sudo ./configure --prefix=$PWD && sudo make && sudo make install

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz && tar -zxvf gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz

COPY v2.1.0.tar.gz ./

RUN tar -zxvf v2.1.0.tar.gz && rm v2.1.0.tar.gz &&\
    mv PASApipeline-2.1.0 PASApipeline && cd PASApipeline && make clean && make && cd .. &&  cp ../conf_files/conf.txt PASApipeline/pasa_conf/ &&\
    cp ../scripts/process_GMAP_alignments_gff3_chimeras_ok.pl PASApipeline/scripts/ &&\
    chmod 775 PASApipeline/scripts/process_GMAP_alignments_gff3_chimeras_ok.pl

COPY augustus.current.tar.gz ./

RUN tar -zxvf augustus.current.tar.gz && rm augustus.current.tar.gz && cd augustus  && make clean && make

COPY Trinity-v2.5.1.tar.gz ./

RUN tar -zxvf Trinity-v2.5.1.tar.gz && \
    mv trinityrnaseq-Trinity-v2.5.1 Trinity &&rm Trinity-v2.5.1.tar.gz && cd Trinity && make && make plugins

RUN git clone https://github.com/alexdobin/STAR.git

COPY v3.0.1.tar.gz ./

RUN tar -zxvf v3.0.1.tar.gz && rm v3.0.1.tar.gz &&\
    cd TransDecoder-3.0.1 && make

COPY gmap-gsnap-2017-11-15.tar.gz ./

RUN tar -zxvf gmap-gsnap-2017-11-15.tar.gz && rm gmap-gsnap-2017-11-15.tar.gz  && \
    mv gmap-2017-11-15  gmap && cd gmap/ && ./configure && make && sudo make install

COPY fasta-36.3.8e.tar.gz ./

RUN tar -zxvf fasta-36.3.8e.tar.gz && rm fasta-36.3.8e.tar.gz &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 && cp /opt/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /usr/local/bin/fasta

COPY BRAKER2.tar.gz ./

RUN tar -zxvf BRAKER2.tar.gz && rm BRAKER2.tar.gz && mv BRAKER* BRAKER

COPY v1.1.1.tar.gz ./

RUN tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i  Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric && sudo cpan -f -i File::Which

RUN mkdir gffread && cd gffread && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /usr/local/bin

COPY genometools-1.5.9.tar.gz ./

RUN tar -zxvf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make

RUN cp ../../code/createUser.py /usr/local/bin && chmod 775 /usr/local/bin/createUser.py

RUN cp ../conf_files/pathToExport.txt /etc/profile.d/pathToExport.sh

RUN cp ../conf_files/extrinsic.M.RM.E.W.P.cfg augustus/config/extrinsic/

RUN rm /opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl

RUN sudo chmod -R 775 /opt/LoReAn/code/

COPY interproscan-5.27-66.0-64-bit.tar.gz ./

RUN tar -pxvzf interproscan-5.27-66.0-64-bit.tar.gz && rm interproscan-5.27-66.0-64-bit.tar.gz

WORKDIR /opt/LoReAn/third_party/software/interproscan-5.27-66.0

RUN mkdir cddblast

COPY ncbi-blast-2.7.1+-x64-linux.tar.gz ./cddblast

RUN cd cddblast && tar -zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz && cp -r ncbi-blast-2.7.1+ ../bin/blast

COPY signalp-4.1f.Linux.tar.gz ./

RUN  tar -xzf signalp-4.1f.Linux.tar.gz -C bin/signalp/4.1 --strip-components 1 && rm signalp-4.1f.Linux.tar.gz

COPY signalp-4.1/signalp bin/signalp/4.1/

RUN mkdir /data_panther

COPY tmhmm-2.0c.Linux.tar.gz ./

RUN  tar -xzf tmhmm-2.0c.Linux.tar.gz -C ./ && cp tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64  bin/tmhmm/2.0c/decodeanhmm && \
     cp tmhmm-2.0c/lib/TMHMM2.0.model  data/tmhmm/2.0c/TMHMM2.0c.model

RUN cp /opt/LoReAn/third_party/conf_files/interproscan.properties ./interproscan.properties

WORKDIR /usr/local/bin

RUN apt-get install -y -q hmmer

COPY trf /usr/bin/

WORKDIR /usr/local

RUN mkdir nseg && cd nseg && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/* && make && mv nseg ../bin && mv nmerge ../bin

COPY RepeatScout-1.0.5.tar.gz ./

RUN tar -xvf RepeatScout* && rm RepeatScout*.tar.gz && mv RepeatScout* RepeatScout && cd RepeatScout && make

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz && \
    tar -xzvf ncbi-rmblastn* && rm ncbi-rmblastn*.tar.gz && mv ncbi-rmblastn*/bin/rmblastn bin && rm -rf ncbi-rmblastn


RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast* && find ncbi-blast* -type f -executable -exec mv {} bin \; &&  rm -rf ncbi-blast*

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i Text::Soundex

COPY RepeatMasker-open-4-0-7.tar.gz ./

RUN tar -xzvf RepeatMasker-open*.tar.gz \
	&& rm -f RepeatMasker-open*.tar.gz && perl -0p -e 's/\/usr\/local\/hmmer/\/usr\/bin/g;' -e 's/\/usr\/local\/rmblast/\/usr\/local\/bin/g;' \
    -e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
    -e 's/TRF_PRGM = ""/TRF_PRGM = "\/usr\/local\/bin\/trf"/g;' RepeatMasker/RepeatMaskerConfig.tmpl > RepeatMasker/RepeatMaskerConfig.pm

RUN cd /usr/local/RepeatMasker && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
	RepeatMasker DateRepeats ProcessRepeats RepeatProteinMask DupMasker util/queryRepeatDatabase.pl \
	util/queryTaxonomyDatabase.pl util/rmOutToGFF3.pl util/rmToUCSCTables.pl

WORKDIR /opt/LoReAn/

RUN chmod a+w /opt/

WORKDIR /data/

CMD /usr/local/bin/createUser.py
