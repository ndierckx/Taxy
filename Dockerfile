FROM debian:latest

ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN apt -y update &&\
    apt -y upgrade &&\
    apt install -y locales &&\
    sed -i /etc/locale.gen -e '/en_[UG].*.UTF-8/s/# //' &&\
    locale-gen &&\
    apt install -y libmce-perl libparallel-forkmanager-perl libdevel-nytprof-perl &&\
    apt install -y mafft ncbi-blast+

COPY ./Taxy0.1.pl /usr/bin

RUN chmod 775 /usr/bin/Taxy0.1.pl

RUN apt clean &&\
    rm -rf /var/lib/apt var/cache/apt &&\
    apt -y purge apt --allow-remove-essential --auto-remove

# podman build -f Dockerfile -t taxy
# podman save -o Taxy.tar taxy
# singularity build --fakeroot Taxy.sif docker-archive://Taxy.tar
# singularity pull docker://ghcr.io/ndierckx/taxy:latest
# ./taxy.sif Taxy0.1.pl
