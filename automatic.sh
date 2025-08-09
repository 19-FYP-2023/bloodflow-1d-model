#!/bin/bash
git checkout local
git pull origin local
# docker build --no-cache -t arteryfe:2017.2.0 .

#loop patient number from 1 to 4374
for i in {10..4374}
#run createCardiacOutput.py  with system aregument for i
do
    
    git pull origin local
    python3 createCardiacOutput.py $i
    if [ $? -ne 0 ]; then
        echo "Error pre processing patient $i"
        exit 1
    fi
    echo "Patient $i processed successfully"
    sleep 1

    git add data/Original_aortic/input.csv
    git commit -m "Processed patient $i"
    git push origin local
    if [ $? -ne 0 ]; then
        echo "Error pushing changes for patient $i"
        exit 1
    fi
    echo "Changes for patient $i pushed successfully"

    docker rm -v -f arteryfe_container
    docker rmi -f arteryfe:2017.2.0
    
    docker build --no-cache -t arteryfe:2017.2.0 .
    docker run --name arteryfe_container -d -it arteryfe:2017.2.0 /bin/bash

    # Wait a moment for container to fully start
    sleep 2

    docker exec arteryfe_container /bin/bash -c "cd bloodflow-1d-model && git checkout local && git pull origin local && chmod +x automaticDocker.sh && ./automaticDocker.sh"
    if [ $? -ne 0 ]; then
        echo "Error running automaticDocker.sh for patient $i"
        exit 1
    fi
    echo "automaticDocker.sh executed successfully for patient $i"

    docker cp arteryfe_container:/home/fenics/bloodflow-1d-model/output/4cycles_last /home/dumindu/modeling/out/output_$i
    if [ $? -ne 0 ]; then
        echo "Error copying output for patient $i"
        exit 1
    fi
    echo "Output copied successfully for patient $i"

    docker stop arteryfe_container

    git checkout skin_model
    git pull origin skin_model

    python3 skin_model_simul.py output_$i

    
done

