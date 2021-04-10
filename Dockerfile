FROM ubuntu:20.10

LABEL author="Clemens Mauksch"

RUN apt-get update && apt-get install sudo

RUN sudo apt-get install -y openjdk-14-jdk default-jre curl python3 cutadapt fastqc trim-galore abyss wget

RUN adduser --disabled-password --gecos '' pipeline
RUN adduser pipeline sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER pipeline
WORKDIR /home/pipeline/
RUN chmod a+rwx /home/pipeline/

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
RUN bash Anaconda3-2020.11-Linux-x86_64.sh -b
RUN rm Anaconda3-2020.11-Linux-x86_64.sh 

ENV PATH /home/pipeline/anaconda3/bin:$PATH

RUN conda update conda
RUN conda update anaconda
RUN conda update --all

RUN conda install -c bioconda bowtie2 samtools