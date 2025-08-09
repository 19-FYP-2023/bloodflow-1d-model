#!/bin/bash

echo "start automaticDocker.sh"


git checkout local
git pull origin local

echo "start setup.py"

sudo python3 setup.py install
if [ $? -ne 0 ]; then
    echo "Error running setup.py"
    exit 1
fi
echo "setup.py executed successfully"

sleep 5
echo "start demo_arterybranch.py"

python3 demo_arterybranch.py config/demo_arterybranch.cfg
if [ $? -ne 0 ]; then
    echo "Error running demo_arterybranch.py"
    exit 1
fi
echo "demo_arterybranch.py executed successfully"

python3 postprocess.py output/4cycles_last/data.cfg