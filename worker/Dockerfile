FROM ubuntu:latest

WORKDIR /proj

RUN apt-get update && apt-get install -y python python-pip octave

ADD worker/requirements.txt /proj

RUN pip install --upgrade pip

RUN pip install -r requirements.txt

ADD python/ python/
ADD testmat/ matlab/
