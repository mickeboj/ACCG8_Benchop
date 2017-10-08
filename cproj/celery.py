from __future__ import absolute_import, unicode_literals
from celery import Celery

app = Celery('cproj',
             broker='amqp://user1:pwd@localhost:5672/vhost1',
             backend='amqp://user1:pwd@localhost:5672/vhost1',
             include=['cproj.tasks'])



if __name__ == '__main__':
    app.start()
