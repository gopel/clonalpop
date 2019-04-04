# Use an official Python runtime as a parent image
FROM python:3.7-slim

# Set the working directory to /ClonalPop
WORKDIR /ClonalPop

# Copy the current directory contents into the container at /ClonalPop
COPY . /ClonalPop

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get -y install gcc
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME World

# Run ClonalPop.py when the container launches
CMD ["python", "ClonalPop.py"]

#Abricate
RUN apt-get update && apt-get -y install sudo \
  emboss \
  bioperl \
  ncbi-blast+ \
  gzip \
  unzip \
  libjson-perl \
  libtext-csv-perl \
  libfile-slurp-perl \
  liblwp-protocol-https-perl \
  libwww-perl \
  git \
  wget && apt-get clean

RUN git clone https://github.com/tseemann/abricate.git
ENV PATH="./abricate/bin:$PATH"
RUN abricate --check
RUN ./abricate/bin/abricate --setupdb
RUN ./abricate/bin/abricate ./abricate/test/assembly.fa

### Breseq  
RUN wget https://github.com/barricklab/breseq/releases/download/v0.33.2/breseq-0.33.2-Linux-x86_64.tar.gz && tar -zxvf breseq-0.33.2-Linux-x86_64.tar.gz && rm -rf v0.8.7.tar.gz
ENV PATH="./breseq-0.33.2-Linux-x86_64/bin:$PATH"

#Spades 
RUN wget https://github.com/ablab/spades/releases/download/v3.13.0/SPAdes-3.13.0-Linux.tar.gz && tar -zxvf SPAdes-3.13.0-Linux.tar.gz && rm -rf SPAdes-3.13.0-Linux.tar.gz
RUN PATH="./SPAdes-3.13.0-Linux/bin:$PATH"

# Linuxbrew
RUN apt-get update \
	&& apt-get install -y --no-install-recommends ca-certificates curl file g++ git locales make uuid-runtime \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& localedef -i en_US -f UTF-8 en_US.UTF-8 \
	&& useradd -m -s /bin/bash linuxbrew \
	&& echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER linuxbrew
WORKDIR /home/linuxbrew
ENV PATH=/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH \
	SHELL=/bin/bash

RUN git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew/Homebrew \
	&& mkdir /home/linuxbrew/.linuxbrew/bin \
	&& ln -s ../Homebrew/bin/brew /home/linuxbrew/.linuxbrew/bin/ \
	&& brew config

# Unicycler
RUN brew install brewsci/bio/unicycler 
RUN brew switch unicycler 0.4.7

# Prokka
RUN brew install brewsci/bio/prokka
RUN brew switch prokka 1.13
