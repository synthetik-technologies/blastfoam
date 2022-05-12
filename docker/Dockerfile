FROM ubuntu:18.04 AS openfoam
RUN apt-get update \
  && apt-get install -y gnupg2 wget software-properties-common \
  && sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -" \
  && add-apt-repository http://dl.openfoam.org/ubuntu \
  && apt-get update
RUN apt-get -y install openfoam9