#!/bin/bash
rabbitmqctl add_user user1 pwd12
rabbitmqctl add_vhost vost1
rabbitmqctl set_permissions -p party-line guest ".*" ".*" ".*"
rabbitmqctl set_user_tags user1 administrator
