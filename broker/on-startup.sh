#!/bin/bash
rabbitmqctl add_user user1 pwd12
rabbitmqctl add_vhost vhost1
rabbitmqctl set_permissions -p vhost1 user1 ".*" ".*" ".*"
rabbitmqctl set_user_tags user1 administrator
