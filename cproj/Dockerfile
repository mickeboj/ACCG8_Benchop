from lordnelsson/testimg

ADD /BENCHOP bench/
ADD /cproj cproj/
ADD celerytest.py ./

RUN apt-get install -y rabbitmq-server

RUN rabbitmqctl add_user user1 pwd
RUN rabbitmqctl add_vhost vhost1
RUN rabbitmqctl set_permissions -p vhost1 user1 ".*" ".*" ".*"
