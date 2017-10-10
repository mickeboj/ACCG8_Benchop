FROM lordnelsson/accbop_base
ADD /BENCHOP bench/
ADD /cproj cproj/
ADD celerytest.py ./
ADD worker/config.py cproj/
CMD celery -A cproj worker -l info
