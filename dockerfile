FROM centos:7

RUN yum -y install epel-release && yum clean all
RUN yum -y install python-pip && yum clean all

# start by building the basic container
MAINTAINER Lauren Thatch <lmthatch@mines.edu>
RUN yum update -y
# add gfortran, debugging tools and make
RUN yum install -y gcc-gfortran gdb make 

RUN yum install -y python-devel

RUN pip install numpy

#COPY primes.f95 /testing/
#COPY pfb_read.f90 /testing/
COPY primes.so /testing/
COPY checkme.py /testing/

WORKDIR /testing/

#CMD f2py -c primes.f95 pfb_read.f90 -m primes
CMD python checkme.py
