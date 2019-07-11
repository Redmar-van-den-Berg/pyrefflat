FROM python:3.6-slim
RUN apt-get update && apt-get install -y git 
RUN git clone https://github.com/sndrtj/pyrefflat
WORKDIR pyrefflat
RUN pip install pytest pyvcf
RUN python setup.py install
