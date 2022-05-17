FROM ubuntu:18.04 AS openfoam
RUN apt-get update \
  && apt-get install -y wget software-properties-common \
  && sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -" \
  && add-apt-repository http://dl.openfoam.org/ubuntu \
  && apt-get update
RUN apt-get -y install openfoam9

FROM ubuntu:18.04 AS build
SHELL [ "/bin/bash", "-lc" ]
COPY --from=openfoam /opt/openfoam9 /opt/openfoam9
RUN apt-get update && apt-get install -y apt-utils build-essential libopenmpi-dev
WORKDIR /opt/blastfoam
COPY . .
RUN echo "source /opt/openfoam9/etc/bashrc" >> ~/.bashrc \
  && echo "source /opt/blastfoam/etc/bashrc" >> ~/.bashrc \
  && source /opt/openfoam9/etc/bashrc \
  && source /opt/blastfoam/etc/bashrc \
  && ./Allwmake -j \
  && ./Allwclean

FROM ubuntu:18.04 as prod
SHELL [ "/bin/bash", "-c" ]
RUN apt-get update && apt-get install -y gnuplot libopenmpi-dev
COPY --from=build /opt /opt
COPY --from=build /root /root
WORKDIR /app/blastfoam
ENTRYPOINT ["/bin/bash", "-ci"]