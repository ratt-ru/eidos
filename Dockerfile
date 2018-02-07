FROM kernsuite/base:3
RUN docker-apt-install python-tk python-pip
ADD . /eidos
RUN pip install /eidos
RUN eidos -h
