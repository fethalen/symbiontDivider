FROM ubuntu:20.10

LABEL author="Clemens Mauksch"

RUN apt-get update --fix-missing && apt-get install sudo
RUN sudo apt-get update --fix-missing

RUN sudo apt-get install -y apt-utils openjdk-14-jdk default-jre curl python3 cutadapt fastqc wget git ncbi-blast+ trim-galore --fix-missing
RUN sudo apt-get install -y pigz

RUN adduser --disabled-password --gecos '' symbiontdivider
RUN adduser symbiontdivider sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER symbiontdivider
WORKDIR /home/symbiontdivider/Docker
RUN sudo  chmod a+rwx /home/symbiontdivider/Docker

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
RUN bash Anaconda3-2020.11-Linux-x86_64.sh -b
RUN rm Anaconda3-2020.11-Linux-x86_64.sh 

ENV PATH /home/symbiontdivider/anaconda3/bin:$PATH

RUN conda update conda
RUN conda update anaconda
RUN conda update --all

RUN conda install -c bioconda bowtie2 megahit

RUN git clone https://github.com/ablab/quast && cd quast && ./setup.py install && cd ..
ENV PATH /home/symbiontdivider/Docker/quast:$PATH
