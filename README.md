# Benchop-As-A-Service

## Applied Cloud Computing Group 8 Project

###Setting up instance
Heat not available for the moment.

Contextualize with CloudInit
Make sure right inputs are for username/password.
Then run contextualize.py(In testing phase right now).

Otherwise use cloud-cfg.txt as configuration in openstack dashboard when creating an instance.


###Installing Docker (in case of crash)
Add GPG key:
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

add Docker rep to apt:
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

Update:
sudo apt-get update

make sure installation is from docker repo:
apt-cache policy docker-ce

Install docker:
sudo apt-get install -y docker-ce

### Docker
In the folder "worker" there is an docker file called dockerfile which is used to contextualize the docker container
### Creating Docker img
$ docker build -t <name> -f <path/to/dockerfile> .

### Starting worker
From /proj directory
$ celery -A cproj worker -l info


### Cleaning scripts docker
sudo docker rm $(sudo docker ps -q -f 'status=exited')
sudo docker rmi $(sudo docker images -q -f "dangling=true")
sudo docker volume rm $(sudo docker volume ls -qf dangling=true)

## Running things
### Starting rabbitmq container
$ sudo docker run --rm --name broker -d  -p 5672:5672 -v ${PWD}/broker:/hooks authentise/rabbitmq

### If no images locally build

#### Flower
sudo docker build -t <name_flower_img> -f flower/Dockerfile .
#### celery worker
sudo docker build -t <name_worker_img> -f worker/Dockerfile .
#### NOTE
Both worker/config.py and flower.config.py must contain appropiate broker address

### Run the containers
sudo docker run -d --rm --name <some_name> -p 5555:5555 <name_flower_img>

sudo docker run -d --rm --name <some_other_name> <name_worker_img>
