#!/bin/bash
rabbitmqctl add_vhost party-line
rabbitmqctl set_permissions -p party-line guest "." "." ".*"
echo "party-line vhost setup"
