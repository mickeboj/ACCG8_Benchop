# http://docs.openstack.org/developer/python-novaclient/ref/v2/servers.html
import time, os, sys
import inspect
from os import environ as env

from  novaclient import client
import keystoneclient.v3.client as ksclient
from keystoneauth1 import loading
from keystoneauth1 import session

flavor = "ssc.small" 
private_net = "SNIC 2017/13-45 Internal IPv4 Network"
floating_ip_pool_name = None
floating_ip = None
image_name = "Ubuntu 16.04 LTS (Xenial Xerus) - latest"

loader = loading.get_plugin_loader('password')

auth = loader.load_from_options(auth_url='https://hpc2n.cloud.snic.se:5000/v3',
                                                                username='s9508',password='****',project_name="SNIC 2017/13-45",project_domain_name='snic',project_id='ad5091c4f42e4defb98eb9550f875f4f',user_domain_name='snic')

sess = session.Session(auth=auth)
nova = client.Client('2.1', session=sess)
print "user authorization completed."

image = nova.glance.find_image(image_name)

flavor = nova.flavors.find(name=flavor)

if private_net != None:
    net = nova.neutron.find_network(private_net)
    nics = [{'net-id': net.id}]
else:
    sys.exit("private-net not defined.")

#print("Path at terminal when executing this file")
#print(os.getcwd() + "\n")
cfg_file_path =  os.getcwd()+'/cloud-cfg.txt'
if os.path.isfile(cfg_file_path):
    userdata = open(cfg_file_path)
else:
    sys.exit("cloud-cfg.txt is not in current working directory")

secgroups = ['default']

print "Creating instance ... "
instance = nova.servers.create(name="vm1", image=image, flavor=flavor, userdata=userdata, nics=nics,security_groups=secgroups)
inst_status = instance.status
print "waiting for 10 seconds.. "
time.sleep(10)

while inst_status == 'BUILD':
    print "Instance: "+instance.name+" is in "+inst_status+" state, sleeping for 5 seconds more..."
    time.sleep(5)
    instance = nova.servers.get(instance.id)
    inst_status = instance.status

print "Instance: "+ instance.name +" is in " + inst_status + "state"
