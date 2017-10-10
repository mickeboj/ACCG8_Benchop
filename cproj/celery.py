from __future__ import absolute_import, unicode_literals
from celery import Celery
from .config import broker

app = Celery('cproj',
             broker=broker,
             backend=broker,
             include=['cproj.tasks'])



if __name__ == '__main__':
    app.start()
