FROM ubuntu:xenial

RUN apt-get clean all && apt-get update && apt-get install -y -q build-essential git wget perl \
    python3.5 python2.7 software-properties-common python3-pip python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2 exonerate \
    bzip2 liblzma-dev libbz2-dev software-properties-common libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev cmake tcsh g++ iputils-ping jellyfish bowtie bioperl apache2 libcairo2-dev libpango1.0-dev libfile-homedir-perl sqlite3 \
    bamtools libbamtools-dev liblzma-dev automake autoconf hmmer

RUN pip install --upgrade pip && pip3 install numpy biopython==1.68 bcbio-gff==0.6.4 pandas==0.19.1 pybedtools==0.7.8 gffutils regex pysam progressbar2 \
    psutil memory_profiler pathlib colorama simplesam tqdm

WORKDIR /opt/

RUN git clone https://github.com/lfaino/LoReAn.git

WORKDIR /opt/LoReAn/third_party/software/

RUN mv trf /usr/local/bin/

RUN tar -zxvf Porechop.tar.gz && cd Porechop && make clean && make && cp porechop/cpp_functions.so  /opt/LoReAn/code/

RUN tar -zxvf AATpackage-r03052011.tgz && rm AATpackage-r03052011.tgz && cd AATpackage-r03052011 && make clean &&\
 sudo ./configure --prefix=$PWD && sudo make && sudo make install

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz && tar -zxvf gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz

RUN tar -zxvf SE-MEI.tar.gz && cd SE-MEI && make

#COPY PASApipeline-v2.3.3.tar.gz ./

RUN wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.3.3/PASApipeline-v2.3.3.tar.gz && \
    tar -zxvf PASApipeline-v2.3.3.tar.gz  && rm PASApipeline-v2.3.3.tar.gz && mv PASApipeline-v2.3.3 PASApipeline && \
    cd PASApipeline && make clean && make && cp ../../scripts/process_GMAP_alignments_gff3_chimeras_ok.pl scripts/ && \
    chmod 775 scripts/process_GMAP_alignments_gff3_chimeras_ok.pl

##    cd .. &&  cp ../conf_files/conf.txt PASApipeline/pasa_conf/ &&\

#RUN apt-get install -y -q bamtools libbamtools-dev liblzma-dev automake autoconf

RUN git clone https://github.com/samtools/htslib.git && cd htslib && autoheader && autoconf && ./configure && make &&\
    sudo make install && cd .. &&  git clone https://github.com/samtools/bcftools.git && cd bcftools && autoheader &&\
  autoconf && ./configure && make && sudo make install && cd .. && git clone https://github.com/samtools/tabix.git &&\
  cd tabix && make && cd .. && git clone https://github.com/samtools/samtools.git && cd samtools && autoheader &&\
   autoconf -Wno-syntax && ./configure && make && sudo make install

RUN apt-get update && apt-get install -y -q --fix-missing libcurl4-gnutls-dev libssl-dev

RUN export TOOLDIR=/opt/LoReAn/third_party/software && git clone https://github.com/Gaius-Augustus/Augustus.git &&\
    mv Augustus augustus && cd augustus  && make clean && make

##COPY Trinity-v2.5.1.tar.gz ./

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.5.tar.gz && tar -zxvf Trinity-v2.8.5.tar.gz && \
    mv trinityrnaseq-Trinity-v2.8.5 Trinity && rm Trinity-v2.8.5.tar.gz && cd Trinity && make && make plugins

RUN git clone https://github.com/lh3/minimap2.git && cd minimap2 && make

RUN git clone https://github.com/alexdobin/STAR.git

##COPY v3.0.1.tar.gz ./

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz &&\
    tar -zxvf salmon-0.14.1_linux_x86_64.tar.gz

RUN wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz && tar -zxvf TransDecoder-v5.5.0.tar.gz && rm TransDecoder-v5.5.0.tar.gz &&\
    mv TransDecoder-TransDecoder-v5.5.0 TransDecoder-v5.5.0 && cd TransDecoder-v5.5.0 && make

##COPY gmap-gsnap-2017-11-15.tar.gz ./

RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz && tar -zxvf gmap-gsnap-2017-11-15.tar.gz && rm gmap-gsnap-2017-11-15.tar.gz  && \
    mv gmap-2017-11-15  gmap && cd gmap/ && ./configure && make && sudo make install

##COPY fasta-36.3.8e.tar.gz ./

RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz && tar -zxvf fasta-36.3.8g.tar.gz && rm fasta-36.3.8g.tar.gz &&\
    cd fasta-36.3.8g/src && make -f ../make/Makefile.linux fasta36 && cp /opt/LoReAn/third_party/software/fasta-36.3.8g/bin/fasta36 /usr/local/bin/fasta

RUN git clone https://github.com/Gaius-Augustus/BRAKER.git && cd BRAKER && chmod -R 777 scripts/ ##&& ln -s /opt/LoReAn/third_party/software/BRAKER/

#COPY v1.1.1.tar.gz ./

RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric && sudo cpan -f -i File::Which && sudo cpan -f -i DBD::SQLite.pm

RUN git clone https://github.com/gpertea/gclib && git clone https://github.com/gpertea/gffread && cd gffread && make release


#git clone https://github.com/gpertea/gclib &&\
#    git clone https://github.com/gpertea/gffread && cd gffread && make && cp ./gffread /usr/local/bin

##COPY genometools-1.5.9.tar.gz ./

RUN wget http://genometools.org/pub/genometools-1.5.9.tar.gz && tar -zxvf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make

RUN cp ../../code/createUser.py /usr/local/bin && chmod 775 /usr/local/bin/createUser.py

RUN rm /opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl

RUN sudo chmod -R 775 /opt/LoReAn/code/

#COPY interproscan-5.27-66.0-64-bit.tar.gz ./

#RUN tar -pxvzf interproscan-5.27-66.0-64-bit.tar.gz && rm interproscan-5.27-66.0-64-bit.tar.gz

#WORKDIR /opt/LoReAn/third_party/software/interproscan-5.27-66.0

#RUN mkdir cddblast

#COPY ncbi-blast-2.7.1+-x64-linux.tar.gz ./cddblast

#RUN cd cddblast && tar -zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz && cp -r ncbi-blast-2.7.1+ ../bin/blast

#COPY signalp-4.1f.Linux.tar.gz ./

#RUN  tar -xzf signalp-4.1f.Linux.tar.gz -C bin/signalp/4.1 --strip-components 1 && rm signalp-4.1f.Linux.tar.gz

#COPY signalp-4.1/signalp bin/signalp/4.1/

#RUN mkdir /data_panther

#COPY tmhmm-2.0c.Linux.tar.gz ./

#RUN  tar -xzf tmhmm-2.0c.Linux.tar.gz -C ./ && cp tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64  bin/tmhmm/2.0c/decodeanhmm && \
#     cp tmhmm-2.0c/lib/TMHMM2.0.model  data/tmhmm/2.0c/TMHMM2.0c.model

#RUN cp /opt/LoReAn/third_party/conf_files/interproscan.properties ./interproscan.properties

WORKDIR /usr/local/bin

#COPY trf ./

WORKDIR /usr/local

RUN mkdir nseg && cd nseg && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/* && make && mv nseg ../bin && mv nmerge ../bin

##COPY RepeatScout-1.0.5.tar.gz ./

RUN wget http://bix.ucsd.edu/repeatscout/RepeatScout-1.0.5.tar.gz && tar -xvf RepeatScout* && rm RepeatScout*.tar.gz && \
    mv RepeatScout* RepeatScout && cd RepeatScout && make

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz && \
    tar -xzvf ncbi-rmblastn* && rm ncbi-rmblastn*.tar.gz && mv ncbi-rmblastn*/bin/rmblastn bin && rm -rf ncbi-rmblastn


RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast* && find ncbi-blast* -type f -executable -exec mv {} bin \; &&  rm -rf ncbi-blast*

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i Text::Soundex

##COPY RepeatMasker-open-4-0-7.tar.gz ./

RUN wget http://www.repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz && tar -xzvf RepeatMasker-open*.tar.gz \
	&& rm -f RepeatMasker-open*.tar.gz && perl -0p -e 's/\/usr\/local\/hmmer/\/usr\/bin/g;' -e 's/\/usr\/local\/rmblast/\/usr\/local\/bin/g;' \
    -e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
    -e 's/TRF_PRGM = ""/TRF_PRGM = "\/usr\/local\/bin\/trf"/g;' RepeatMasker/RepeatMaskerConfig.tmpl > RepeatMasker/RepeatMaskerConfig.pm

RUN cd /usr/local/RepeatMasker && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
	RepeatMasker DateRepeats ProcessRepeats RepeatProteinMask DupMasker util/queryRepeatDatabase.pl \
	util/queryTaxonomyDatabase.pl util/rmOutToGFF3.pl util/rmToUCSCTables.pl

RUN chmod -R 777 RepeatMasker/

WORKDIR /opt/LoReAn/

RUN cp /opt/LoReAn/code/lorean /usr/local/bin && chmod 775 /usr/local/bin/lorean

RUN apt-get install -y locales && locale-gen en_US.UTF-8  && update-locale

RUN chmod a+w /opt/

WORKDIR /data/

ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/opt/LoReAn/third_party/software/augustus/bin:/opt/LoReAn/third_party/software/LoReAn/:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1:/opt/LoReAn/third_party/software/gm_et_linux_64/gmes_petap:/opt/LoReAn/third_party/software/iAssembler-v1.3.2.x64:/opt/LoReAn/third_party/software/PASApipeline/scripts:/opt/LoReAn/third_party/software/PASApipeline/:/opt/LoReAn/third_party/software/Trinity:/opt/LoReAn/third_party/software/AATpackage-r03052011/bin:/opt/LoReAn/third_party/software/LoReAn/third_party:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc:/opt/LoReAn/third_party/scripts:/opt/LoReAn/code/:/opt/LoReAn/third_party/software/STAR/bin/Linux_x86_64:/opt/LoReAn/third_party/software/PASApipeline/bin:/opt/LoReAn/third_party/software/genometools-1.5.9/bin:/opt/LoReAn/third_party/software/BRAKER/scripts/:/opt/LoReAn/third_party/software/interproscan-5.27-66.0:/opt/LoReAn/third_party/software/TransDecoder-3.0.1:/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util:/opt/LoReAn/third_party/conf_files:/opt/:/usr/local/RepeatMasker:/usr/local/RepeatScout:/opt/LoReAn/third_party/software/PASApipeline/misc_utilities/:/opt/LoReAn/third_party/software/minimap2/:/opt/LoReAn/third_party/software/SE-MEI/:/opt/LoReAn/third_party/software/salmon-latest_linux_x86_64/bin/:/opt/LoReAn/third_party/software/gffread"

ENV AUGUSTUS_CONFIG_PATH="/opt/LoReAn/third_party/software/augustus/config/"

ENV GENEMARK_PATH="/opt/LoReAn/third_party/software/gm_et_linux_64/gmes_petap/"

ENV PERL5LIB="/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/PerlLib/:/opt/LoReAn/third_party/scripts/:/opt/LoReAn/third_party/software/PASApipeline/PerlLib/"

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/SW/src/"

ENV SAMTOOLS_PATH="/usr/local/bin/"

ENV BLAST_PATH="/usr/local/bin/"

ENV TMPDIR="/tmp"


CMD /usr/local/bin/createUser.py
