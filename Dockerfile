FROM ubuntu:20.04 AS openfoam
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y wget software-properties-common \
  && sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -" \
  && add-apt-repository http://dl.openfoam.org/ubuntu \
  && apt-get update
RUN apt-get -y install openfoam12

FROM ubuntu:20.04 as prod
SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive
COPY --from=openfoam /opt/openfoam12 /opt/openfoam12
COPY blastfoam_*.deb .
RUN apt-get -qq -o=Dpkg::Use-Pty=0 update && apt-get -qq -o=Dpkg::Use-Pty=0 install -y gnuplot libopenmpi-dev
RUN dpkg -i blastfoam_*.deb
RUN echo "source /opt/openfoam12/etc/bashrc" >> ~/.bashrc \
  && echo "source /opt/blastfoam/etc/bashrc" >> ~/.bashrc \
  && source ~/.bashrc
WORKDIR /app/blastfoam
ENTRYPOINT ["/bin/bash", "-ci"]
