# Benchop-As-A-Service

## Applied Cloud Computing Group 8 Project

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
