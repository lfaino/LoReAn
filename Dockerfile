FROM ubuntu:bionic

RUN apt-get clean all && apt-get update && apt-get install -y -q --fix-missing build-essential git wget perl \
    python3.6 python2.7 software-properties-common python3-pip python-pip debconf-utils sudo python-numpy cmake samtools bedtools zlib1g-dev libc6 aptitude \
    libdbd-mysql-perl libdbi-perl libboost-all-dev libncurses5-dev bowtie default-jre parallel nano bowtie2 exonerate \
    bzip2 liblzma-dev libbz2-dev software-properties-common libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev cmake tcsh g++ iputils-ping jellyfish bowtie bioperl apache2 libcairo2-dev libpango1.0-dev libfile-homedir-perl sqlite3 \
    bamtools libbamtools-dev liblzma-dev automake autoconf libssl-dev libmysqlclient-dev mysql-client libsqlite3-dev libmysql++-dev \
    libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libboost-iostreams-dev zlib1g-dev libbamtools-dev libbz2-dev \
    liblzma-dev libncurses5-dev libssl-dev libcurl3-dev libboost-all-dev hmmer

RUN pip install --upgrade pip && pip3 install numpy==1.17.1

RUN pip3 install biopython==1.73 bcbio-gff==0.6.4 pandas==0.19.1 \
    pybedtools==0.7.8 gffutils==0.9 regex==2019.8.19 pysam==0.15.3 progressbar2==3.43.1 \
    psutil==5.6.3 memory_profiler==0.55.0 pathlib==1.0.1 colorama==0.4.1 simplesam==0.1.3 tqdm==4.35.0 \
    argcomplete==1.10.0 argh==0.26.2 ordereddict==1.1 pycurl==7.43.0 pyfaidx==0.5.5.2 pygobject python-apt \
    python-dateutil==2.8.0 python-utils==2.3.0 pytz==2019.2 simplejson==3.16.0 six==1.12.0 unattended-upgrades==0.1

WORKDIR /opt/

RUN git clone https://github.com/lfaino/LoReAn.git

WORKDIR /opt/LoReAn/third_party/software/

RUN wget --no-check-certificate http://github.com/bbuchfink/diamond/releases/download/v0.9.26/diamond-linux64.tar.gz && tar -zxvf diamond-linux64.tar.gz

RUN mv trf /usr/local/bin/

RUN tar -zxvf Porechop.tar.gz && cd Porechop && make clean && make && cp porechop/cpp_functions.so  /opt/LoReAn/code/

RUN tar -zxvf iAssembler-v1.3.2.x64.tgz && rm iAssembler-v1.3.2.x64.tgz && tar -zxvf gm_et_linux_64_4.48_3.60.tar.gz && rm gm_et_linux_64_4.48_3.60.tar.gz

RUN git clone --recursive https://github.com/dpryan79/SE-MEI.git && cd SE-MEI && make

RUN wget --no-check-certificate https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.3.3/PASApipeline-v2.3.3.tar.gz && \
    tar -zxvf PASApipeline-v2.3.3.tar.gz  && rm PASApipeline-v2.3.3.tar.gz && mv PASApipeline-v2.3.3 PASApipeline && \
    cd PASApipeline && make clean && make && cp ../../scripts/process_GMAP_alignments_gff3_chimeras_ok.pl scripts/ && \
    chmod 775 scripts/process_GMAP_alignments_gff3_chimeras_ok.pl

RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && bunzip2 htslib-1.9.tar.bz2 && tar -xvf htslib-1.9.tar && \
    mv htslib-1.9 htslib && cd htslib && autoheader && autoconf && ./configure && make &&\
    sudo make install && cd ..

RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && bunzip2 bcftools-1.9.tar.bz2 && \
    tar -xvf bcftools-1.9.tar && mv bcftools-1.9 bcftools && cd bcftools && autoheader &&\
    autoconf && ./configure && make && sudo make install && cd ..

RUN git clone https://github.com/samtools/tabix.git && cd tabix && make && cd ..

RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && bunzip2 samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar && mv samtools-1.9 samtools && cd samtools && autoheader && \
    autoconf -Wno-syntax && ./configure && make && sudo make install

RUN export TOOLDIR=/opt/LoReAn/third_party/software && git clone https://github.com/Gaius-Augustus/Augustus.git &&\
    mv Augustus augustus && cd augustus  && make clean && make

RUN wget --no-check-certificate https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.5.tar.gz && tar -zxvf Trinity-v2.8.5.tar.gz && \
    mv trinityrnaseq-Trinity-v2.8.5 Trinity && rm Trinity-v2.8.5.tar.gz && cd Trinity && make && make plugins

RUN wget --no-check-certificate https://github.com/lh3/minimap2/archive/v2.17.tar.gz && tar -zxvf v2.17.tar.gz && mv minimap2-2.17 minimap2 && cd minimap2 && make

RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz && tar -zxvf 2.7.3a.tar.gz && mv STAR-2.7.3a STAR

RUN wget --no-check-certificate https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz &&\
    tar -zxvf salmon-0.14.1_linux_x86_64.tar.gz

RUN wget --no-check-certificate https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz && tar -zxvf TransDecoder-v5.5.0.tar.gz && rm TransDecoder-v5.5.0.tar.gz &&\
    mv TransDecoder-TransDecoder-v5.5.0 TransDecoder-v5.5.0 && cd TransDecoder-v5.5.0 && make

RUN wget --no-check-certificate http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz && tar -zxvf gmap-gsnap-2017-11-15.tar.gz && rm gmap-gsnap-2017-11-15.tar.gz  && \
    mv gmap-2017-11-15  gmap && cd gmap/ && ./configure && make && sudo make install

RUN wget --no-check-certificate http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz && tar -zxvf fasta-36.3.8g.tar.gz && rm fasta-36.3.8g.tar.gz &&\
    cd fasta-36.3.8g/src && make -f ../make/Makefile.linux fasta36 && cp /opt/LoReAn/third_party/software/fasta-36.3.8g/bin/fasta36 /usr/local/bin/fasta

RUN wget --no-check-certificate  https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.5.tar.gz && tar -zxvf v2.1.5.tar.gz &&\
    mv BRAKER-2.1.5 BRAKER && cd BRAKER && chmod -R 777 scripts/ ##&& ln -s /opt/LoReAn/third_party/software/BRAKER/


RUN wget --no-check-certificate https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar -zxvf v1.1.1.tar.gz && rm v1.1.1.tar.gz

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i YAML && sudo cpan -f -i Hash::Merge && sudo cpan -f -i  Logger::Simple && sudo cpan -f -i Parallel::ForkManager &&\
    sudo cpan -f -i Config::Std && sudo cpan -f -i Scalar::Util::Numeric && sudo cpan -f -i File::Which && sudo cpan -f -i DBD::SQLite.pm

RUN git clone https://github.com/gpertea/gclib && git clone https://github.com/gpertea/gffread && cd gffread && make release

RUN wget --no-check-certificate http://genometools.org/pub/binary_distributions/gt-1.6.0-Linux_x86_64_x86_64-64bit.tar.gz \
    && tar -zxvf gt-1.6.0-Linux_x86_64_x86_64-64bit.tar.gz && mv gt-1.6.0-Linux_x86_64_x86_64-64bit genometools

RUN cp ../../code/createUser.py /usr/local/bin && chmod 775 /usr/local/bin/createUser.py

RUN rm /opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl

RUN sudo chmod -R 775 /opt/LoReAn/code/

RUN wget --no-check-certificate http://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz && tar -zxvf gth-1.7.3-Linux_x86_64-64bit.tar.gz

RUN git clone  https://github.com/gpertea/cdbfasta.git && cd cdbfasta/ && make

RUN rm *.tar.gz

WORKDIR /usr/local/bin

WORKDIR /usr/local

RUN mkdir nseg && cd nseg && wget --no-check-certificate ftp://ftp.ncbi.nih.gov/pub/seg/nseg/* && make && mv nseg ../bin && mv nmerge ../bin

RUN wget --no-check-certificate http://bix.ucsd.edu/repeatscout/RepeatScout-1.0.5.tar.gz && tar -xvf RepeatScout* && rm RepeatScout*.tar.gz && \
    mv RepeatScout* RepeatScout && cd RepeatScout && make

RUN wget --no-check-certificate ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz && \
    tar -xzvf ncbi-rmblastn* && rm ncbi-rmblastn*.tar.gz && mv ncbi-rmblastn*/bin/rmblastn bin && rm -rf ncbi-rmblastn


RUN wget --no-check-certificate ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast* && find ncbi-blast* -type f -executable -exec mv {} bin \; &&  rm -rf ncbi-blast*

RUN sudo perl -MCPAN -e shell && sudo cpan -f -i Text::Soundex

RUN wget --no-check-certificate http://www.repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz && tar -xzvf RepeatMasker-open*.tar.gz \
	&& rm -f RepeatMasker-open*.tar.gz && perl -0p -e 's/\/usr\/local\/hmmer/\/usr\/bin/g;' -e 's/\/usr\/local\/rmblast/\/usr\/local\/bin/g;' \
    -e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
    -e 's/TRF_PRGM = ""/TRF_PRGM = "\/usr\/local\/bin\/trf"/g;' RepeatMasker/RepeatMaskerConfig.tmpl > RepeatMasker/RepeatMaskerConfig.pm

RUN cd /usr/local/RepeatMasker && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
	RepeatMasker DateRepeats ProcessRepeats RepeatProteinMask DupMasker util/queryRepeatDatabase.pl \
	util/queryTaxonomyDatabase.pl util/rmOutToGFF3.pl util/rmToUCSCTables.pl

RUN chmod -R 777 RepeatMasker/

RUN cp /opt/LoReAn/third_party/software/diamond /usr/local/bin/ && chmod 777 /usr/local/bin/diamond

WORKDIR /opt/LoReAn/

RUN cp /opt/LoReAn/code/lorean /usr/local/bin/ && chmod -R 775 /usr/local/bin/ && chmod a+w /opt/

RUN apt-get install -y locales && locale-gen en_US.UTF-8  && update-locale

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fastqToFa && mv fastqToFa /usr/local/bin && chmod 775 /usr/local/bin/fastqToFa

WORKDIR /data/

RUN mkdir /home/lorean

ENV HOME="/home/lorean"

ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/opt/LoReAn/third_party/software/:/opt/LoReAn/third_party/software/augustus/bin:/opt/LoReAn/third_party/software/LoReAn/:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1:/opt/LoReAn/third_party/software/gm_et_linux_64:/opt/LoReAn/third_party/software/iAssembler-v1.3.2.x64:/opt/LoReAn/third_party/software/PASApipeline/scripts:/opt/LoReAn/third_party/software/PASApipeline/:/opt/LoReAn/third_party/software/Trinity:/opt/LoReAn/third_party/software/AATpackage-r03052011/bin:/opt/LoReAn/third_party/software/LoReAn/third_party:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc:/opt/LoReAn/third_party/scripts:/opt/LoReAn/code/:/opt/LoReAn/third_party/software/STAR/bin/Linux_x86_64:/opt/LoReAn/third_party/software/PASApipeline/bin:/opt/LoReAn/third_party/software/genometools/bin:/opt/LoReAn/third_party/software/BRAKER/scripts/:/opt/LoReAn/third_party/software/interproscan-5.27-66.0:/opt/LoReAn/third_party/software/TransDecoder-3.0.1:/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util:/opt/LoReAn/third_party/conf_files:/opt/:/usr/local/RepeatMasker:/usr/local/RepeatScout:/opt/LoReAn/third_party/software/PASApipeline/misc_utilities/:/opt/LoReAn/third_party/software/minimap2/:/opt/LoReAn/third_party/software/SE-MEI/:/opt/LoReAn/third_party/software/salmon-latest_linux_x86_64/bin/:/opt/LoReAn/third_party/software/gffread:/opt/LoReAn/third_party/software/gth-1.7.3-Linux_x86_64-64bit/bin:/opt/LoReAn/third_party/software/cdbfasta/:/opt/RepeatMasker:/opt/RepeatMasker/util:/opt/RepeatModeler:/opt/RepeatModeler/util:/opt/coseg"

ENV AUGUSTUS_BIN_PATH="/opt/LoReAn/third_party/software/augustus/bin/"

ENV AUGUSTUS_SCRIPTS_PATH="/opt/LoReAn/third_party/software/augustus/scripts/"

ENV DIAMOND_PATH="/usr/local/bin/"

ENV CDBTOOLS_PATH="/opt/LoReAn/third_party/software/cdbfasta/"

ENV AUGUSTUS_CONFIG_PATH="/opt/LoReAn/third_party/software/augustus/config/"

ENV GENEMARK_PATH="/opt/LoReAn/third_party/software/gm_et_linux_64"

ENV PERL5LIB="/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/PerlLib/:/opt/LoReAn/third_party/scripts/:/opt/LoReAn/third_party/software/PASApipeline/PerlLib/"

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/SW/src/"

ENV SAMTOOLS_PATH="/usr/local/bin/"

ENV BLAST_PATH="/usr/local/bin/"

ENV TMPDIR="/tmp"

ENV LANG=C

CMD /usr/local/bin/lorean
