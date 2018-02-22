FROM kernsuite/base:3
RUN docker-apt-install python-pip
ADD . /eidos
RUN pip install -U pip pyyaml
RUN pip install /eidos
RUN eidos -h
