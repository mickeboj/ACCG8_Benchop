FROM lordnelsson/accbop_base

RUN pip install flask

ADD /BENCHOP bench/
ADD /cproj cproj/
ADD flask/app.py ./
ADD worker/config.py cproj/

EXPOSE 5000

CMD python app.py
