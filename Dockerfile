FROM ubuntu:xenial

WORKDIR /opt/

RUN apt-get clean all && apt-get update && apt-get install -y -q --fix-missing build-essential git wget perl \
    python3.5 python2.7 software-properties-common python3-pip python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2 exonerate \
    bzip2 liblzma-dev libbz2-dev software-properties-common libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev cmake tcsh g++ iputils-ping jellyfish bowtie bioperl apache2 libcairo2-dev libpango1.0-dev \
    bamtools libbamtools-dev liblzma-dev automake autoconf \
    libcurl4-gnutls-dev libssl-dev \
    hmmer locales && \
    locale-gen en_US.UTF-8  && update-locale && \
    echo "mysql-server mysql-server/root_password password lorean" | debconf-set-selections && \
    echo "mysql-server mysql-server/root_password_again password lorean" | debconf-set-selections && \
    apt-get install -y mysql-server mysql-client mysql-common && \
    pip3 install numpy biopython==1.68 bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils regex pysam matplotlib progressbar2 \
    psutil memory_profiler pathlib colorama simplesam tqdm Flask && \
    chmod a+w /opt/ && \
    git clone https://github.com/lfaino/LoReAn.git

WORKDIR /opt/LoReAn/third_party/software

COPY PASApipeline-v2.3.3.tar.gz ./
COPY Trinity-v2.5.1.tar.gz ./
COPY v3.0.1.tar.gz ./
COPY gmap-gsnap-2017-11-15.tar.gz ./
COPY fasta-36.3.8e.tar.gz ./
COPY v1.1.1.tar.gz ./
COPY genometools-1.5.9.tar.gz ./
COPY interproscan-5.27-66.0-64-bit.tar.gz ./
COPY ncbi-blast-2.7.1+-x64-linux.tar.gz ./
COPY signalp-4.1f.Linux.tar.gz ./
COPY tmhmm-2.0c.Linux.tar.gz ./
COPY trf /usr/local/bin
COPY RepeatScout-1.0.5.tar.gz /usr/local
COPY RepeatMasker-open-4-0-7.tar.gz /usr/local

RUN tar -zxvf Porechop.tar.gz && cd Porechop && make clean && make && cp porechop/cpp_functions.so /opt/LoReAn/code/ && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean && \
    sudo ./configure --prefix=$PWD && sudo make && sudo make install && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz && tar -zxvf gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz && \
    tar -zxvf SE-MEI.tar.gz && cd SE-MEI && make && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf PASApipeline-v2.3.3.tar.gz  && rm PASApipeline-v2.3.3.tar.gz &&\
    mv PASApipeline-v2.3.3 PASApipeline && cd PASApipeline && make clean && make && cd .. &&  cp ../conf_files/conf.txt PASApipeline/pasa_conf/ &&\
    cp ../scripts/process_GMAP_alignments_gff3_chimeras_ok.pl PASApipeline/scripts/ && chmod 775 PASApipeline/scripts/process_GMAP_alignments_gff3_chimeras_ok.pl && \
    cd /opt/LoReAn/third_party/software && \
    git clone https://github.com/samtools/htslib.git && cd htslib && autoheader && autoconf && ./configure && make &&\
    sudo make install && cd .. &&  git clone https://github.com/samtools/bcftools.git && cd bcftools && autoheader &&\
    autoconf && ./configure && make && sudo make install && cd .. && git clone https://github.com/samtools/tabix.git &&\
    cd tabix && make && cd .. && git clone https://github.com/samtools/samtools.git && cd samtools && autoheader &&\
    autoconf -Wno-syntax && ./configure && make && sudo make install && \
    cd /opt/LoReAn/third_party/software && \
    export TOOLDIR=/opt/LoReAn/third_party/software && git clone https://github.com/Gaius-Augustus/Augustus.git &&\
    mv Augustus augustus && cd augustus  && make clean && make && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf Trinity-v2.5.1.tar.gz && \
    mv trinityrnaseq-Trinity-v2.5.1 Trinity &&rm Trinity-v2.5.1.tar.gz && cd Trinity && make && make plugins && \
    cd /opt/LoReAn/third_party/software && \
    git clone https://github.com/lh3/minimap2.git && cd minimap2 && make && \
    cd /opt/LoReAn/third_party/software && \
    git clone https://github.com/alexdobin/STAR.git && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf v3.0.1.tar.gz && rm v3.0.1.tar.gz &&\
    cd TransDecoder-3.0.1 && make && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf gmap-gsnap-2017-11-15.tar.gz && rm gmap-gsnap-2017-11-15.tar.gz  && \
    mv gmap-2017-11-15  gmap && cd gmap/ && ./configure && make && sudo make install && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf fasta-36.3.8e.tar.gz && rm fasta-36.3.8e.tar.gz &&\
    cd fasta-36.3.8e/src && make -f ../make/Makefile.linux fasta36 && cp /opt/LoReAn/third_party/software/fasta-36.3.8e/bin/fasta36 /usr/local/bin/fasta && \
    cd /opt/LoReAn/third_party/software && \
    git clone https://github.com/Gaius-Augustus/BRAKER.git && \
    tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz && \
    sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i  Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric && sudo cpan -f -i File::Which && \
    mkdir gffread && cd gffread && git clone https://github.com/gpertea/gclib &&\
    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /usr/local/bin && \
    cd /opt/LoReAn/third_party/software && \
    tar -zxvf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make && \
    cd /opt/LoReAn/third_party/software && \
    cp ../../code/createUser.py /usr/local/bin && chmod 775 /usr/local/bin/createUser.py && \
    rm /opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl && \
    sudo chmod -R 775 /opt/LoReAn/code/ && \
    tar -pxvzf interproscan-5.27-66.0-64-bit.tar.gz && rm interproscan-5.27-66.0-64-bit.tar.gz && \
    cd /opt/LoReAn/third_party/software/interproscan-5.27-66.0 && \
    mkdir cddblast && \
    mv /opt/LoReAn/third_party/software/ncbi-blast-2.7.1+-x64-linux.tar.gz cddblast/ && \
    cd cddblast && tar -zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz && cp -r ncbi-blast-2.7.1+ ../bin/blast && \
    cd /opt/LoReAn/third_party/software/interproscan-5.27-66.0 && \
    mv /opt/LoReAn/third_party/software/signalp-4.1f.Linux.tar.gz ./ && \
    tar -xzf signalp-4.1f.Linux.tar.gz -C bin/signalp/4.1 --strip-components 1 && rm signalp-4.1f.Linux.tar.gz && \
    cp signalp-4.1/signalp bin/signalp/4.1/ && \
    mkdir /data_panther && \
    cd /opt/LoReAn/third_party/software/interproscan-5.27-66.0 && \
    mv /opt/LoReAn/third_party/software/tmhmm-2.0c.Linux.tar.gz ./ && \
    tar -xzf tmhmm-2.0c.Linux.tar.gz -C ./ && cp tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64  bin/tmhmm/2.0c/decodeanhmm && \
    cp tmhmm-2.0c/lib/TMHMM2.0.model  data/tmhmm/2.0c/TMHMM2.0c.model && \
    cp /opt/LoReAn/third_party/conf_files/interproscan.properties ./interproscan.properties

WORKDIR /usr/local

RUN mkdir nseg && cd nseg && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/* && make && mv nseg ../bin && mv nmerge ../bin && \
    cd /usr/local && \
    tar -xvf RepeatScout* && rm RepeatScout*.tar.gz && mv RepeatScout* RepeatScout && cd RepeatScout && make && \
    cd /usr/local && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz && \
    tar -xzvf ncbi-rmblastn* && rm ncbi-rmblastn*.tar.gz && mv ncbi-rmblastn*/bin/rmblastn bin && rm -rf ncbi-rmblastn && \
    cd /usr/local && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast* && find ncbi-blast* -type f -executable -exec mv {} bin \; &&  rm -rf ncbi-blast* && \
    sudo perl -MCPAN -e shell && sudo cpan -f -i Text::Soundex && \
    cd /usr/local && \
    tar -xzvf RepeatMasker-open*.tar.gz && \
	  rm -f RepeatMasker-open*.tar.gz && perl -0p -e 's/\/usr\/local\/hmmer/\/usr\/bin/g;' -e 's/\/usr\/local\/rmblast/\/usr\/local\/bin/g;' \
    -e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
    -e 's/TRF_PRGM = ""/TRF_PRGM = "\/usr\/local\/bin\/trf"/g;' RepeatMasker/RepeatMaskerConfig.tmpl > RepeatMasker/RepeatMaskerConfig.pm && \
    cd /usr/local/RepeatMasker && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
  	RepeatMasker DateRepeats ProcessRepeats RepeatProteinMask DupMasker util/queryRepeatDatabase.pl \
  	util/queryTaxonomyDatabase.pl util/rmOutToGFF3.pl util/rmToUCSCTables.pl && \
    chmod -R 777 RepeatMasker/

WORKDIR /data/

CMD /usr/local/bin/createUser.py
