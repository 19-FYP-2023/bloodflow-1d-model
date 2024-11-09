FROM quay.io/fenicsproject/stable:current

USER root

RUN sudo apt-get update
RUN sudo apt-get -y install git python3-setuptools
RUN sudo apt-get autoremove

# WORKDIR /bloodflow
# RUN mkdir /bloodflow

# COPY . .
# COPY ./arteryfe /bloodflow/arteryfe



# RUN chmod +x ./app.py

# RUN python3 app.py
RUN git clone https://github.com/19-FYP-2023/bloodflow-1d-model.git

RUN rm -f /home/fenics/bloodflow-1d-model/data/example_inlet.csv
COPY Input/example_inlet.csv /home/fenics/bloodflow-1d-model/data

RUN cd bloodflow-1d-model && sudo python3 setup.py install && cd ..

# RUN cd bloodflow-1d-model && sudo python3 demo_arterybranch.py config/demo_arterybranch.cfg && cd ..

# RUN cd bloodflow-1d-model && sudo python3 demo_arterybranch.py config/demo_arterybranch.cfg && cd ..

# RUN cd bloodflow-1d-model && sudo python3 postprocess.py output/4cycles_last/data.cfg && cd ..

# RUN cd bloodflow-1d-model && sudo python3 postprocess.py output/4cycles_last/data.cfg && cd ..

# RUN cd bloodflow-1d-model && python3 demo_arterybranch.py --cfg config/demo_arterybranch.cfg && cd ..

# CMD ["python3","/bloodflow/app.py"]




# RUN cd bloodflow 
# RUN sudo python3 setup.py install 
# RUN cd ..