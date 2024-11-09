FROM dolfinx/dolfinx:latest

USER root

RUN apt-get update && apt-get install -y \
    git \
    python3-setuptools \
    && apt-get autoremove -y \
    && apt-get clean -y

# Clone the repository
RUN git clone https://github.com/19-FYP-2023/bloodflow-1d-model.git

# Remove the existing example_inlet.csv and copy the new one
RUN rm -f /home/fenics/bloodflow-1d-model/data/example_inlet.csv
COPY Input/example_inlet.csv /home/fenics/bloodflow-1d-model/data

# Install the Python package
RUN cd bloodflow-1d-model && python3 setup.py install && cd ..

# Uncomment and modify the following lines if you need to run specific scripts
# RUN cd bloodflow-1d-model && python3 demo_arterybranch.py config/demo_arterybranch.cfg && cd ..
# RUN cd bloodflow-1d-model && python3 postprocess.py output/4cycles_last/data.cfg && cd ..

# Set the working directory
WORKDIR /home/fenics/bloodflow-1d-model
