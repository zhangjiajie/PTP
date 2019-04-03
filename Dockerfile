FROM python:3.5-slim
MAINTAINER Jiajie Zhang <bestzhangjiajie@gmail.com>

COPY . /app
RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir -r /app/requirements.txt

WORKDIR /app
RUN python3 setup.py install
