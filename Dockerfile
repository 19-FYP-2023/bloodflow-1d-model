FROM quay.io/fenicsproject/stable:2017.2.0

USER root

RUN sudo apt-get update
RUN sudo apt-get -y install git python3-setuptools
RUN sudo apt-get autoremove

WORKDIR /bloodflow
# RUN mkdir /bloodflow

COPY . .
# COPY ./arteryfe /bloodflow/arteryfe



# RUN chmod +x ./app.py

# RUN python3 app.py

CMD ["python3","/bloodflow/app.py"]



# RUN git clone https://github.com/akdiem/bloodflow.git
# RUN cd bloodflow 
# RUN sudo python3 setup.py install 
# RUN cd ..
