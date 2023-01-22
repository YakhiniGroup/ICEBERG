FROM python:3.8-slim as builder

COPY ./requirements.txt /root/

WORKDIR /root/


RUN apt update \
&& apt install -y libmariadb-dev \
gcc \
libcogl-pango-dev \
libcairo2-dev \
libtool \
linux-headers-amd64 \
musl-dev \
python-jinja2 \
libffi-dev \
libssl-dev \
libjpeg-dev \
zlib1g-dev \
python3-dev \
python3-pip \
python3-numpy \
# libncurses6 libncursesw6 \
python3-pandas

RUN apt-get update && apt-get install --no-install-recommends -y \
libncurses5-dev \
libbz2-dev \
liblzma-dev \
libcurl4-gnutls-dev \
zlib1g-dev \
libssl-dev \
wget \
make \
perl \
bzip2 \
gnuplot \
ca-certificates \
gawk && \
apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN apt update && apt install libncurses6

RUN apt-get update && apt-get install -y bwa
RUN apt-get update && apt-get install -y samtools

RUN pip install --upgrade pip && \
pip install --no-cache -r requirements.txt


FROM python:3.8-slim
COPY --from=builder /usr /usr
COPY --from=builder /lib /lib
COPY --from=builder /lib64 /lib64

COPY iceberg /root/iceberg
COPY images/iceberg-logo.jpg /root/images/iceberg-logo.jpg

# COPY DOCKER-TEST/TIGIT/1 /root/EXPERIMENTS_FOLDERS/TIGIT/1
# COPY DOCKER-TEST/MOCK /root/EXPERIMENTS_FOLDERS/MOCK
# COPY DOCKER-TEST/GENOME /root/GENOME
# COPY DOCKER-TEST/TIGIT_1_INPUT.yaml /root/TIGIT_1_INPUT.yaml



# CMD ['python','iceberg/analyzer']
# command to run on container start
# CMD [ "python", "-m", "unittest", "./test/test_iceberg.py" ]
# CMD [ "sleep", '10' ]

# CMD ["tail", "-f", "/dev/null"]